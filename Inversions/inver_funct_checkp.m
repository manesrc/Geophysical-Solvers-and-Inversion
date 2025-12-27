function [inverted_variable,error_vect,phi_t,accepted_samples,rejected_samples,time_inversion,relevDOF, InfoMesh, InfoProblem] = ...
    inver_funct_checkp(frwdmethod, inversion_case, cond_dir, Long_X, maxDepth, nel_x, nel_y, sigma_d, len_prop, LAB_definition, nsteps,SynthTemperature, ...
    plot_cut_elements,velo01,LABpriorGeometry, make_up, file_name_save)

    %% 0. RESUME CHECK
    checkpoint_file = file_name_save;
    resume_flag = false;
    
    if exist(checkpoint_file, 'file')
        fprintf('Resuming from checkpoint: %s\n', checkpoint_file);
        load(checkpoint_file);
        resume_flag = true;
    else
        fprintf('Starting new inversion: %s\n', checkpoint_file);
        i_start = 2; 
    end

    %% Path Management
    if ~contains(path, '00_general_operations')
        % Add other folders if the primary one isn't found
        addpath('00_general_operations', '01_plots', '02_mesh_cond_SF', ...
            '03_LevelSet', '04_TempInMesh', '05_Stokes_2D', '89_LAB_definition', '07_inversion_definitions');
        if frwdmethod == 2
            addpath '06_Forward2ndMethod'
        else
            addpath '06_Forward1stMethod'
        end 
    end

    if ~resume_flag
        % inputs         
        InfoProblem.L_ref = max(Long_X,maxDepth);
        [InfoProblem, InfoMaterial] = Inputs_file(InfoProblem);
        InfoMesh = generate_mesh_related_data(nel_x, nel_y,Long_X,maxDepth,InfoProblem);
        
        % compute prior distribution
        plot_prior_idea = 0; % set to 1 to plot prior distribution 
        [prior_bounds_phi, relevDOF] = generate_prior(LABpriorGeometry.LAB_minDepth, ...
            LABpriorGeometry.LAB_maxDepth, LABpriorGeometry.LAB_XChanges, InfoMesh.X, maxDepth, plot_prior_idea);
        
        % First LAB proposal m^0
        curve_case = LAB_definition.case;
        varargin.rescaleY = LAB_definition.Y_0;
        varargin.data_name = 'datos_isotherm';
        varargin.model_width = Long_X; 
        varargin.model_height = maxDepth; 
        varargin.y0 = LAB_definition.depth_ini;
        varargin.mx = LAB_definition.mx ; 
        InfoLAB0 = defineLAB(curve_case,varargin);
        % estimate parameters for Temp field
        grad_T1 = (InfoProblem.T_LAB-InfoProblem.T_sup)/(maxDepth - mean(InfoLAB0.LABy)); 
        InfoProblem.k2 = InfoProblem.k1 * grad_T1 / InfoProblem.grad_apprT2; % estimate conductivity in Asthenosphere
        InfoProblem.q2 = InfoProblem.k2 * InfoProblem.grad_apprT2; % estimate flux in Asthenosphere
        InfoProblem.T_inf = InfoProblem.T_LAB + InfoProblem.grad_apprT2 * mean(InfoLAB0.LABy);
        
        %% Compute temperature and pressure estimation to compute the matrix
        varargin.plotLS = 0; varargin.nel_x = InfoMesh.nel_x; varargin.nel_y = InfoMesh.nel_y;
        InfoMesh.LS_new = LevelSet(InfoMesh.X,InfoLAB0,varargin);
        if make_up
            InfoMesh.LS_new = LSreinit(InfoMesh.X,InfoMesh.T,InfoMesh.LS_new')';
        end 
        h_caract = InfoProblem.L_ref * max(abs(InfoMesh.fin_x-InfoMesh.ini_x)/InfoMesh.nel_x,abs(InfoMesh.fin_y-InfoMesh.ini_y)/InfoMesh.nel_y);
        InfoProblem.tolerance = h_caract * 0.02;
        points2add2LAB = 1000; 
        [InfoLAB0.LABx,InfoLAB0.LABy] = increase_data_LAB(InfoLAB0,points2add2LAB);
        [InfoMesh.list_Omega1,InfoMesh.list_Omega2,InfoMesh.list_cut, InfoMesh.list_edge1,InfoMesh.list_edge2, interphasePoints] = findElemNitsche_LS(InfoProblem.tolerance, ...
            plot_cut_elements,InfoMesh,InfoLAB0,InfoProblem);
        
        [Temp_est,pres_est] = estimate_Temp_pres(InfoMesh.X,InfoMesh.nel_x,InfoMesh.nel_y,InfoMesh.LS_new, InfoMaterial,InfoProblem);
        
        %% compute problem results for m^0
        plot_velo.fig = 0;
        plot_velo.LAB = [interphasePoints(1,:); interphasePoints(2,:)] / InfoProblem.L_ref; 
        plot_velo.parameters = 0; 
        
        % compute initial temperature for m^0
        if frwdmethod == 2
            plot_velo.fig_velo = 0; 
            % first iteration solves Stokes problem:
            [u_st0, ~] = ComputeStokesProblem(Temp_est,pres_est,InfoMesh,InfoProblem,InfoMaterial,plot_velo);
            nonlin = 1; 
            reduction = 17; % choose reduction level for basis
            name_basis = ['98_data_base/newUbasis',num2str(InfoMesh.nel_x),'_',num2str(InfoMesh.nel_y),'.mat'];
            basis1 = load(name_basis,'U_basis');
            if ((reduction == 17) || (reduction == 10))
                [U_basis,velo_comp] = sum_square_velocities(basis1.U_basis,InfoMesh.nel_x,InfoMesh.nel_y);
            else
                U_basis = basis1.U_basis;
            end
            [~,SigmaEigVal,~] = svd(U_basis);
            flag_lin_dep = sum(diag(SigmaEigVal)/SigmaEigVal(1,1) < 1e-8);
            U_trimmed = update_basis_data(U_basis,velo_comp,InfoMesh.list_Omega1,[]);
            alpha_old = []; 
            flag_basis_changes = [];
            [Temp1, u_mant1, alpha1] = method2forwardsolver(Temp_est, pres_est, alpha_old, flag_basis_changes, u_st0, nonlin, velo_comp, U_trimmed, flag_lin_dep, cond_dir,InfoMesh,InfoProblem,InfoMaterial);
        else
            varargin.vel_01 = velo01;
            [Temp1, results1,InfoMesh1,InfoProblem1] = method1forwardsolver(Temp_est, pres_est,InfoMesh,InfoProblem,InfoMaterial,varargin, plot_velo);
        end 
        
        
        % define inversion case and data
        matrixSupFlux = []; % Initialize empty
        dof2compare = [];
        
        switch inversion_case
            case 'Temp_Omega'
                direct_results = Temp1;
                Data4inversion = SynthTemperature(:,1); 
                variance = 100;
            case 'Temp_sup'
                dof2compare = find(InfoMesh.X(:,2)>=0.75*max(InfoMesh.X(:,2)));
                Temperature_results_Omega = Temp1;
                direct_results = Temperature_results_Omega(dof2compare);
                Data4inversion = SynthTemperature(dof2compare,1); 
                variance = 100;
            case 'flux_sup'
                elements_considered = (size(InfoMesh.T,1) - InfoMesh.nel_x +1):1:size(InfoMesh.T,1); 
                Temperature_results_Omega = Temp1;
                matrixSupFlux = generateMatrixFluxSup(InfoProblem.k1, elements_considered,InfoMesh);
                direct_results = matrixSupFlux * Temperature_results_Omega;
                Data4inversion = matrixSupFlux * SynthTemperature(:,1); 
                variance = 1e-8; 
            otherwise
                error('unknown inversion type')
        end 
        
        
        error0 = (1/(2*variance)) * norm(direct_results - Data4inversion)^2; 
        phi_t = InfoMesh.LS_new; 
        error_vect = zeros(nsteps,1);
        inverted_variable = zeros(length(relevDOF),nsteps);
        error_vect(1) = error0;
        inverted_variable(:,1) = phi_t(relevDOF); 
        
        % save initial velocities and optimized U_basis combination
        if frwdmethod == 2
            u_st_acc = u_st0;
            alpha_it_acc = alpha1;
            Omega1_old_acc = InfoMesh.list_Omega1; 
        end 


        %% DEFINE INVERSE SET
        accepted_samples = 0;
        rejected_samples = 0;
        prior_rejection = 0;
        
    end % End of Resume Flag check

    %% Main Loop
    t1 = tic;
    save_interval = max(1, floor(nsteps / 10)); 
    
    for i = i_start:nsteps
        
        % --- Periodic Checkpoint ---
        if mod(i, save_interval) == 0
            fprintf('Saving iteration %d/%d...\n', i, nsteps);
            i_start = i; 
            save(checkpoint_file); % Save workspace
        end

        % elements to draw from
        list_elem2draw = [InfoMesh.list_cut; InfoMesh.list_edge1];
        randomPick_LScomp2change = randi(length(list_elem2draw(:,1)));
        elem_sample = list_elem2draw(randomPick_LScomp2change,1);
        Te = InfoMesh.T(elem_sample,:);
        Xe = InfoMesh.X(Te,:); 
        c1 = Xe(1,:) + (Xe(3,:) - Xe(1,:)) / 2; 
        delta_h = sigma_d * randn(1);
        phi_e_t = phi_t(Te);
        delta_phi = computeDeltaPhi(delta_h, c1,len_prop,phi_e_t, Xe,InfoMesh);
        phi_int = phi_t + delta_phi;
        if make_up 
            phi_int = LSreinit(InfoMesh.X, InfoMesh.T, phi_int')';
        end 
        passess_prior = (phi_int >= prior_bounds_phi(:,1) & phi_int <= prior_bounds_phi(:,2)); 
        
        if sum(passess_prior) == length(phi_int)
            InfoMesh.LS_new = phi_int;
            [InfoMesh.list_Omega1,InfoMesh.list_Omega2,InfoMesh.list_cut, InfoMesh.list_edge1,InfoMesh.list_edge2, interphasePoints] = findElemNitsche_LS(InfoProblem.tolerance,plot_cut_elements,InfoMesh,[],InfoProblem);
            plot_velo.LAB = [interphasePoints(1,:); interphasePoints(2,:)] / InfoProblem.L_ref; 
            assert(size(InfoMesh.list_edge1,1) == size(InfoMesh.list_edge2,1))
            [Temp_est,pres_est] = estimate_Temp_pres(InfoMesh.X,InfoMesh.nel_x,InfoMesh.nel_y,InfoMesh.LS_new, InfoMaterial,InfoProblem);
            grad_T1 = (InfoProblem.T_LAB-InfoProblem.T_sup)/mean(interphasePoints(2,:)); 
            InfoProblem.k2 = InfoProblem.k1 * grad_T1 / InfoProblem.grad_apprT2;      
            InfoProblem.q2 = InfoProblem.k2 * InfoProblem.grad_apprT2; 
            InfoProblem.T_inf = InfoProblem.T_LAB + InfoProblem.grad_apprT2 * (maxDepth - mean(interphasePoints(2,:)) );
            
            if frwdmethod == 2
                [U_basis_new, flag1] = update_basis_data(U_basis,velo_comp,InfoMesh.list_Omega1,Omega1_old_acc);
                [Temp_it, u_st_it, alpha_it] = method2forwardsolver(Temp_est, pres_est, alpha_it_acc, flag1, u_st_acc, nonlin, velo_comp, U_basis_new, flag_lin_dep, cond_dir,InfoMesh,InfoProblem,InfoMaterial);
            else
                [Temp_it, ~,~,~] = method1forwardsolver(Temp_est, pres_est, InfoMesh,InfoProblem,InfoMaterial,varargin,plot_velo);
            end 
            
            % Use switch-case optimized logic (pre-computed matrices)
            if strcmp(inversion_case, 'Temp_Omega')
                direct_results_it = Temp_it;
            elseif strcmp(inversion_case, 'Temp_sup')
                direct_results_it = Temp_it(dof2compare);
            elseif strcmp(inversion_case, 'flux_sup')
                direct_results_it = matrixSupFlux * Temp_it;
            end
            
            error_int = (1/(2*variance)) * norm(direct_results_it - Data4inversion)^2; 
            u_random = rand();
            logu_rnd = log(u_random);
            ratio_log = -error_int + error_vect(i-1);
            if ratio_log >= logu_rnd
                inverted_variable(:,i) = phi_int(relevDOF); 
                error_vect(i) = error_int; 
                accepted_samples = accepted_samples+ 1; 
                phi_t = phi_int;
                if frwdmethod == 2 
                    alpha_it_acc = alpha_it;
                    u_st_acc = u_st_it;
                    Omega1_old_acc = InfoMesh.list_Omega1; 
                end 
            else
                inverted_variable(:,i) = inverted_variable(:,i-1);
                error_vect(i) = error_vect(i-1); 
                rejected_samples = rejected_samples + 1; 
            end 
        else
            inverted_variable(:,i) = inverted_variable(:,i-1);
            error_vect(i) = error_vect(i-1); 
            prior_rejection = prior_rejection + 1;
        end
    end 
    time_inversion = toc(t1);
    
    % Clean up checkpoint file after successful completion 
    delete(checkpoint_file); 
end


function MatrixComputingFluxSup = generateMatrixFluxSup(k,elements,InfoMesh)
    % this function takes the vector of temperatures in the mesh and calculates
    % the fluxes for those elements in "elements" list
    % inputs:
    % temp : vector of nodal temperatures in the mesh
    % elements : relevant elements to compute the flux 
    % InfoMesh : mesh information
    X = InfoMesh.X;
    T = InfoMesh.T; 
    num_dof = size(X,1); 
    num_elem = length(elements);
    num_dim = size(X,2); 
    num_rows = num_elem * num_dim;
    % compute fluxes
    MatrixComputingFluxSup = zeros(num_rows,num_dof); % 1 per elements x 2 dimensions; num_dof Temp solution
    % compute the fluxes at the elements center
    center_isopar_coord = [0 0];
    [~,Nxi,Neta] = shapeFunctions(InfoMesh.elemType, size(T,2), center_isopar_coord);
    for ii = 1:num_elem
        elem_ii = elements(ii);
        Te = T(elem_ii,:);
        Xe = X(Te,:);
        % compute derivative of shape functions in X,Y
        jacob = [Nxi*Xe(:,1)  Nxi*Xe(:,2);
                    Neta*Xe(:,1) Neta*Xe(:,2)]; 
        dNxy = jacob\[Nxi; Neta]; 
        ind_1 = 2*(ii-1) +1;
        ind_2 = ind_1 + 1; 
        MatrixComputingFluxSup(ind_1,Te) = -k * dNxy(1,:);
        MatrixComputingFluxSup(ind_2,Te) = -k * dNxy(2,:);
    end 
end 