function [T_combined, results,InfoMesh,InfoProblem] = method1forwardsolver(Temp_est, pres_est, InfoMesh,InfoProblem,InfoMaterial,varargin2,plot_up)
    % some information needed to runthis code:
    InfoMesh.LS_mesh = InfoMesh.LS_new; 
    InfoMesh.list1 = InfoMesh.list_Omega1;
    InfoMesh.list2 = InfoMesh.list_Omega2;
    InfoMesh.Delaunay = 1; 
    InfoMesh.ngp_interface = 9; 
    InfoProblem.c =2; 
    T_ref = 1;
    InfoMesh.elemType_LAB = 0; 
    InfoProblem.q2 = InfoProblem.k2*InfoProblem.grad_apprT2;
    InfoProblem.vel_01 = varargin2.vel_01;
    InfoMaterial.T_ref0 = 273;
    InfoMesh.ngp_Nits = 16;
    
    %% gauss points and shape functions (save structure) (Section inputs)
    [MeshIntegration.pospg,MeshIntegration.pespg] = quadrature(InfoMesh.elemType,InfoMesh.nne); 
    [MeshIntegration.N,MeshIntegration.Nxi,MeshIntegration.Neta]  = shapeFunctions(InfoMesh.elemType,InfoMesh.nne, MeshIntegration.pospg); 
    
    %% define dof in sub-domains
    [DOF1,DOF2,DOF_interphase,DOF_G] = define_dof(InfoMesh);
    
    try
        % First, check if DOF_interphase contains only valid indices
        %disp('Checking DOF_interphase for invalid entries...');
        if any(DOF_interphase <= 0) || any(mod(DOF_interphase, 1) ~= 0)
            DOF_interphase
            %mod(DOF_interphase,1)
            error('DOF_interphase contains invalid indices (non-positive or non-integers).');
        end
        % Now proceed with the original calculation
        DOF_int1 = DOF_interphase(InfoMesh.LS_mesh(DOF_interphase) >= 0);
        DOF_int2 = DOF_interphase(InfoMesh.LS_mesh(DOF_interphase) < 0);
    catch ME
        disp('An error occurred:');
        disp(ME.message);
    
        % Check the variables causing the issue
        disp('Problematic variables:');
        disp('DOF_interphase:');
        disp(DOF_interphase);
        
        % Try to extract InfoMesh.LS_mesh for these indices to see what the error is
        try
            disp('LS(DOF_interphase):');
            disp(InfoMesh.LS_mesh(DOF_interphase));  % This might fail
        catch innerME
            disp('Error while trying to access InfoMesh.LS_mesh(DOF_interphase):');
            disp(innerME.message);
        end
    
        % Check InfoMesh.LS_mesh itself
        disp('InfoMesh.LS_mesh:');
        disp(InfoMesh.LS_mesh);
    
        % Rethrow the original error to see the complete error stack
        rethrow(ME);
    end
    
    
    %% STEP 1: Solution in Omega_1
    problem_int = 1; 
    % Obtain FEM matrices
    addpath '06_Forward1stMethod/04_BulkAndNitscheMatrices'
    [K_fem1,f_fem1,Ar_phys_1] = K_FEM_n_sparse(InfoMesh,MeshIntegration, InfoProblem,problem_int);
    % Nitsche matrix computation
    [G1,M1,b1,m1,~] = Nitsche_matrices2D_sparse(InfoMesh, InfoProblem,problem_int,Ar_phys_1);  
    N1 = -G1 + M1;
    vn1 = -b1+m1;
    % FEM + Nitsche
    MAT1 = K_fem1+N1;
    vect1 = f_fem1 + vn1;
    % solve OMEGA_1
    %t1 = tic;
    T_1_nd = MAT1(DOF1,DOF1)\vect1(DOF1);
    %toc(t1)
    
    % dimensional solution OMEGA_1 [in K]
    T_1_dim = T_ref*T_1_nd; 
    
    %% STEP 2: Solution in Omega_2 considering mean flux on Gamma Bottom
    problem_int = 2;
    % Obtain FEM matrices
    [K_fem2,f_fem2,Ar_phys_2] = K_FEM_n_sparse(InfoMesh,MeshIntegration, InfoProblem,problem_int);
    % Nitsche matrix computation
    [G2,M2,b2,m2,~] = Nitsche_matrices2D_sparse(InfoMesh, InfoProblem,problem_int,Ar_phys_2);  
    [g2_mean_short,line_int_vect] = compute_g_mean(DOF_G,InfoMesh,InfoProblem); 
    g2_mean = [g2_mean_short; zeros(size(K_fem2,1)-length(g2_mean_short),1)];
    N2 = -G2 + M2;
    vn2 = -b2+m2;
    % FEM + Nitsche
    MAT2 = K_fem2+N2;
    vect2_gmean = f_fem2 + g2_mean + vn2;
    % solve OMEGA_2
    %t2 = tic;
    T_2_nd_gmean = MAT2(DOF2,DOF2)\vect2_gmean(DOF2);
    %toc(t2)
    % % dimensional solution OMEGA_1 [in K]
    T_2_dim_gmean = T_ref*T_2_nd_gmean; 
    
    %% Stokes problem
    %plot_up.LAB = [InfoLAB.LABx InfoLAB.LABy];
    
    addpath '06_Forward1stMethod/06_stokes_2D_split'
    Temp1 = zeros(size(InfoMesh.X,1),1);
    Temp2 = zeros(size(InfoMesh.X,1),1);
    %warning('for the inverse problem we need to use the temperature estimation and not invert matrix twice')
    %Temp1(DOF1)  = Temp_est(DOF1);
    %Temp2(DOF2) = Temp_est(DOF2);
    Temp1(DOF1)  = T_1_dim;
    Temp2(DOF2) = T_2_dim_gmean;
    %plot_up1.fig = 1;
    %plot_up1.parameters = 1;
    
    %varargin2.u_st = []; varargin2.vel_01= 1;
    
    if InfoProblem.vel_01 == 1
        %t_st = tic;
        [u_st,p_st] = ComputeStokesProblem(Temp1,Temp2,InfoMesh,InfoProblem,InfoMaterial,plot_up); % in [m/s] and [Pa]
        %toc(t_st)
    else
        u_st = zeros(size(InfoMesh.X_v));
        p_st = zeros(size(InfoMesh.X_v,1),1);
    end 
    
    % DISCLAIMER: u_st is originally quadratic, but it's converted to linear before taking it out to use afterwards
    
    %% Advection-diffusion T_2 in Omega_2
    %u_nd = u_st * ((L_ref*InfoMaterial.rho_ref)/(T_ref*k_ref))^(1/3);
    %InfoMaterial.calorific = InfoMaterial.calorific_dim * ((L_ref^2 * T_ref * InfoMaterial.rho_ref^2)/ (k_ref^2))^(1/3);
    
    addpath '06_Forward1stMethod/07_convection_gradT'
    
    if InfoProblem.vel_01 == 1
        % convection matrix:
        G2_matrix = computeConvectionMatrix(u_st, Temp_est, pres_est,MeshIntegration,InfoMesh,InfoProblem,InfoMaterial);
        MAT2_new = K_fem2(DOF2,DOF2)+N2(DOF2,DOF2)+G2_matrix(DOF2,DOF2);
        vect2_indep = f_fem2(DOF2) + vn2(DOF2) + g2_mean(DOF2) ;
        T2_2_nd_gmean_adv = MAT2_new\vect2_indep;
    else
        MAT2_new = MAT2(DOF2,DOF2);
        vect2_indep = f_fem2(DOF2) + vn2(DOF2) + g2_mean(DOF2) ;
        T2_2_nd_gmean_adv = T_2_nd_gmean;
    end 
    
    
    
    %% Compute gradient matrices
    % compute gradient matrix such that flux_1 = k_1 * \nabla T_1 = G_1 * T_1
    
    % num_points_per_elem = 5;
    
    addpath '06_Forward1stMethod/08_Minim_Operations'
    figurePoints = 0;       % plot where the points end 
    % problem_int = 1;    % which matrix to compute
    % min_ratio = 0.2;     % tolerance of are inside the physical domain note as both domains are 
    %                               % considered the complementary ratio comp_ratio = 0.8 is also considered                              
    [PosMat1,PosMatrix1,matG1,PosMatrix2,PosMatrix21,matG2] = compute_Grad_interphase2(DOF1, DOF2, figurePoints, InfoProblem, InfoMesh);
    InfoMesh.PosMat = PosMat1;
    % problem_int = 2;
    % [matG2,~] = compute_Grad_interphase2(num_points_per_elem,DOF2,problem_int,figurePoints,InfoProblem,InfoMesh,Ar_phys_2,min_ratio);
    %% just to compare
    addpath '06_Forward1stMethod/10_rest_minim'
    %                                                           (MAT2_new2,vect2_indep,g2_mean,line_int_vect,matG1,matG2,T_1_nd,T_2_gmean,DOF2,DOF_G,InfoProblem,InfoMesh)
    results = function_comparisson_rest(MAT2_new,vect2_indep,g2_mean(DOF_G),line_int_vect,matG1,matG2,T_1_nd,T2_2_nd_gmean_adv,DOF2,DOF_G,InfoProblem,InfoMesh);
    Temp2(DOF2) = results.T2nd_res;
    
    results.g_mean = g2_mean;
    results.T1nd = T_1_nd;
    results.u_mantle = u_st;
    results.p_mantle = p_st;
    results.DOF1 = DOF1;
    results.DOF2 = DOF2;
    results.DOF_interphase = DOF_interphase;
    results.DOF_G = DOF_G;
    results.PosMat = PosMat1;
    
    results.DOF_int1 = DOF_int1; 
    results.DOF_int2 = DOF_int2; 
    Temp1(DOF_int2) = 0; 
    Temp2(DOF_int1) = 0; 
    T_combined = Temp1+Temp2; 
    results.T_combined = T_combined;
end 