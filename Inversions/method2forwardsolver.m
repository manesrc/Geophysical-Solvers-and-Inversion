function [Temp_min_vect, u_min, alpha_min] = method2forwardsolver(Temp_est, pres_est, alpha_old, flag1, u_old, nonlin, velo_comp, U_basis, ...
    flag_dep, cond_dir, InfoMesh, InfoProblem, InfoMaterial)
% Solve the forward problem for temperature and velocity
% INPUTS:
% - Temp_est, pres_est: Estimated temperature and pressure
% - alpha_old: Previous iteration coefficients for the velocity basis
% - flag1: Logical flag for basis components (two columns: new and old flags)
% - u_old: Previous velocity field
% - U_basis: Current velocity basis
% - flag_dep: Flag to determine dependency on SVD
% OUTPUTS:
% - Temp_min_vect: Minimization result for temperature
% - u_min: Velocity field solution
% - alpha_min: Minimization result for velocity coefficients

    % Define the basis to use
    if flag_dep == 0
        % Least-squares solution for alpha
        alpha_LS = (U_basis' * U_basis) \ (U_basis' * u_old(:)); 
        U_ast = U_basis; 
    else
        % Reduced basis using SVD
        [U, S, ~] = svd(U_basis, 'econ');
        svd_threshold = 1e-8;
        n_col_SVD = sum(diag(S) / S(1, 1) > svd_threshold);
        U_ast = U(:, 1:n_col_SVD);
        alpha_LS = U_ast' * u_old(:);
    end 
    
    % Initialize or resize alpha_old2
    if isempty(alpha_old)
        % Initialize with alpha_LS if no previous data
        alpha_old2 = alpha_LS; 
    else
        % Map old alpha coefficients to the new basis using flag1
        % flag1 = [components2trim_new components2trim_old]
        indices1 = 1:size(flag1, 1);
        % Indices of components to keep in old and new alpha
        indices_in_old_alpha = indices1(flag1(:, 2) == 0); % Old alpha
        indices_in_new_alpha = indices1(flag1(:, 1) == 0); % New alpha
        % generate this variable in full size
        fake_alpha_old = zeros(size(flag1,1),1);
        % first initialize all values with the Least-Sq approx
        fake_alpha_old(indices_in_new_alpha) = alpha_LS; 
        % for those that there's information replace with the previously obtained solution
        fake_alpha_old(indices_in_old_alpha) = alpha_old;
        % recover the valuable component values 
        alpha_old2 = fake_alpha_old(indices_in_new_alpha);
    end 

   % t1 = tic;
    [u_min, alpha_min, ~, ~, ~, ~, ~] = run_matlab_minim_vect_new_velo(alpha_LS, alpha_old2, nonlin,...
        cond_dir, U_ast, InfoMesh, InfoProblem, InfoMaterial,Temp_est,pres_est);
    %toc(t1)
    Temp_min_vect = solveTemperatureProblem(u_min,Temp_est,pres_est,cond_dir,InfoMesh,InfoProblem, ...
        InfoMaterial);
end 

