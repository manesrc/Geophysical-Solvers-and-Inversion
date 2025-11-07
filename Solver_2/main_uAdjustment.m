% Goal: Given a target Lithosphere-Astosphere-Boundary (LAB) depth, find a
% temperature field and convection velocities that are consistent with it.
% This script solves a coupled conduction-convection temperature
% equation and a Stokes equation for fluid flow.
clear; clc
%% Model Properties
% ----------------------------------------------------------------------
% NOTE TO USER:
% Pre-computed velocity bases exist for specific model setups to speed up execution.
%
% 1. For 'Afonso08' case:
%    Use: Long_X = 3000*1000; maxDepth = 400*1000;
%         nel_x = 80; nel_y = 12;
%    (You will also need to uncomment the 'check if there exists results'
%    block below to load the reference 'Afonso08' temperature/velocity)
%
% 2. For 'linear', 'sinusoidal', or 'matrix_data' cases:
%    Use: Long_X = 660*1000; maxDepth = 660*1000;
%         nel_x = 30; nel_y = 30;
%    (And ensure str_name_basis is set to '98_data_base/OldUbasis30_30.mat')
% ----------------------------------------------------------------------
Long_X = 3000*1000; % [m] X-Direction (Model Width)
maxDepth = 400*1000; % [m] Y-Direction (Model Depth)
InfoProblem.L_ref = max(Long_X,maxDepth); % [m] Reference length for non-dimensionalization
% Add path for helper functions
addpath '00_general_operations'
% Load physical constants, material properties, and problem parameters
[InfoProblem, InfoMaterial] = Inputs_file(InfoProblem);
%% Generate Mesh (Fixed): domain for background meshes
InfoMesh.ini_x = 0; InfoMesh.fin_x = Long_X / InfoProblem.L_ref; 
InfoMesh.ini_y = 0; InfoMesh.fin_y = maxDepth/InfoProblem.L_ref; 
% Elements for Temperature/Pressure (Q1)
InfoMesh.nel_x = 80;        % Number of elements in x-direction
InfoMesh.nel_y = 12;        % Number of elements in y-direction
InfoMesh.elemType = 1;      % 1 = linear quadrilateral element
InfoMesh.nne = 4;           % Number of nodes per element
addpath '02_mesh_cond_SF'
% Create T/P mesh: X = Nodal Coords, T = Element Connectivity
[InfoMesh.X,InfoMesh.T] = CreaMalla_rectangulo(InfoMesh.nne, InfoProblem.L_ref, InfoMesh);
% Elements for Velocity (Q2)
InfoMesh.ngp_pois=9;        % Number of Gauss points (for Stokes)
InfoMesh.nne_v = 9;         % Number of nodes per element (9-node quadratic)
[InfoMesh.X_v,InfoMesh.T_v] = CreaMalla_rectangulo(InfoMesh.nne_v, InfoProblem.L_ref, InfoMesh);
%% Define initial LAB distribution
addpath '89_LAB_definition'
% --- Inputs for defineLAB function ---
% 'curve_case' specifies which LAB geometry to use.
% 'varargin' passes a structure with parameters needed for that case.
curve_case = 'Afonso08';    % Options: 'linear', 'sinusoidal', 'matrix_data', 'Afonso08'
varargin.model_width = Long_X;
varargin.model_height = maxDepth; 
% varargin.data_name = '98_data_base/datos_isotherm.mat'; % Used by 'matrix_data' case
% varargin.mx = -0.1;                        % Used by 'linear' case
% varargin.y0 =  0.9*maxDepth;               % Used by 'linear' case
% --------------------------------------
InfoLAB = defineLAB(curve_case,varargin);
% Calculate initial thermal gradients based on the new LAB depth
grad_T1 = (InfoProblem.T_LAB-InfoProblem.T_sup)/(maxDepth - mean(InfoLAB.LABy)); % [K/m]
InfoProblem.k2 = InfoProblem.k1 * grad_T1 / InfoProblem.grad_apprT2; % [W/(m K)]        
disp(['conductivity k_2 = ',num2str(round(InfoProblem.k2,2))])
% Set thermal boundary conditions
cond_dir = 1; % 1 = Dirichlet (fixed T), 0 = Neumann (fixed flux)
if cond_dir == 1
    InfoProblem.T_inf = InfoProblem.T_LAB + InfoProblem.grad_apprT2 * mean(InfoLAB.LABy);
else
    InfoProblem.q2 = InfoProblem.k2 * InfoProblem.grad_apprT2; % [W/(m^2)]
end 
%% check if there exists results or execute the convergence code
% This block is for the 'Afonso08' case, which has a reference solution
if strcmp(curve_case, 'Afonso08') == 1
    str_name_convergence = ['98_data_base/Afonso08_conv_nelx_',num2str(InfoMesh.nel_x),'_',num2str(InfoMesh.nel_y),'_TLAB_1573K.mat'];
    if exist(str_name_convergence,'file') == 2
        load(str_name_convergence);
    else
        addpath '06_OptimizationProcedure'
        % generate a converged state
        [InfoLAB,InfoMesh,InfoProblem,InfoMaterial,u_mant,p_omega,Temp] = generate_converged_state(cond_dir,InfoLAB,InfoMesh,InfoMaterial,InfoProblem);
        save(str_name_convergence,'InfoLAB','InfoMesh','InfoProblem','InfoMaterial','u_mant','p_omega','Temp')
    end 
    str_name_basis = ['98_data_base/Ubasis',num2str(InfoMesh.nel_x),'_',num2str(InfoMesh.nel_y),'.mat'];
else
    str_name_basis = ['98_data_base/OldUbasis',num2str(InfoMesh.nel_x),'_',num2str(InfoMesh.nel_y),'.mat'];
end
%% velocity basis: check or compute
% The velocity basis is a set of pre-computed flow fields (modes) used for model order reduction.
addpath '03_LevelSet'
% Define path to the pre-computed basis file
if exist(str_name_basis,'file') == 2
    disp('Loading existing velocity basis...');
    load(str_name_basis)
else
    disp('Velocity basis not found. Computing new one...');
    % define plots to check:
    addpath '05_Stokes_2D'
    plot_settings.plotLS = 0; plot_settings.parameters = 0;
    plot_settings.fig_velo = 1; plot_settings.nel_x = InfoMesh.nel_x; 
    plot_settings.nel_y = InfoMesh.nel_y; plot_settings.tolPlot = InfoProblem.L_ref*0.01; 
    
    % Define material properties *specifically for the basis computation*
    InfoMaterial2.rho_depth_temp = 0; 
    InfoMaterial2.rho1 = 3300; InfoMaterial2.rho2 = 3000; % [kg/m^3]
    InfoMaterial2.mu1 = 1e21; InfoMaterial2.mu2 = 1e21; % [Pa*s]
    InfoMaterial2.rho_ref = max(InfoMaterial2.rho1,InfoMaterial2.rho2);
    InfoMaterial2.mu_ref = max(InfoMaterial2.mu1,InfoMaterial2.mu2);
    InfoMaterial2.rho_pres_temp = 0; InfoMaterial2.rho_discont = 1; %,1
    InfoMaterial2.mu_discont = 1; InfoMaterial2.mu_pres_temp = 0; %,1
    InfoMaterial2.mu_smooth = 0; InfoMaterial.rho_alphaT = 0; 
    InfoMaterial2.gravity_units = InfoMaterial.gravity_units;
    % compute basis
    viscosity_decayWDist = 1; % (1 yes, 0 no)
    U_basis = VelocityBasis(str_name_basis, plot_settings,InfoProblem,InfoMaterial2,InfoMesh,viscosity_decayWDist);
end
%% Compute temperature and pressure estimation for matrix computation
varargin.plotLS = 0; varargin.nel_x = InfoMesh.nel_x; varargin.nel_y = InfoMesh.nel_y;
% Create Level Set function (0 at LAB, <0 above, >0 below)
InfoMesh.LS_val = LevelSet(InfoMesh.X,InfoLAB,varargin);
InfoMesh.LS_new = InfoMesh.LS_val;
% Estimate T and P with a simple bi-linear distribution (initial guess)
[Temp_est,pres_est] = estimate_Temp_pres(InfoMesh.X,InfoMesh.nel_x,InfoMesh.nel_y,InfoMesh.LS_val, InfoMaterial,InfoProblem);
% Find elements cut by the LAB interface for Nitsche's method
InfoProblem.tol1 = 0.01*max(Long_X/InfoMesh.nel_x,maxDepth/InfoMesh.nel_y);
plot1 = 1; 
[InfoMesh.list_Omega1, InfoMesh.list_Omega2, InfoMesh.list_cut,InfoMesh.list_edge1,InfoMesh.list_edge2]= findElemCrossedLS(InfoProblem.tol1,plot1,InfoMesh,InfoLAB,InfoProblem);
%% Compute reference velocity field
% Since a converged state was not loaded, we solve the Stokes problem
% using the estimated T/P fields to get a reference velocity 'u_mant'.
do_plots_check.plot = 1; 
if do_plots_check.plot == 1
    do_plots_check.LAB(:,1) = InfoLAB.LABx/InfoProblem.L_ref;  
    do_plots_check.LAB(:,2) = InfoLAB.LABy/InfoProblem.L_ref; 
    do_plots_check.mesh = InfoMesh.X_v;
end 
do_plots_check.fig_velo = 1; 
addpath '05_Stokes_2D'
do_plots_check.parameters = 0; 
[u_mant, p1] = ComputeStokesProblem(Temp_est,pres_est,InfoMesh,InfoProblem,InfoMaterial,do_plots_check);
%% SVD over velocity basis
% Perform Singular Value Decomposition on the basis to get a reduced set of principal 
% components (modes) for velocity.
analyze_basis = 1;
modify_basis.execute = 0; 
modify_basis.n2sum = []; 
% 'reduction' flag controls how the basis is truncated/modified
reduction = 17; % 17 = "from base grouped by four take out those from Lithosphere"
addpath '100_Ubasis_changes'
U_basis_norm = normalize_Ubasis(U_basis);
do_plots_check.plot = 1;
do_plots_check.LABx = InfoLAB.LABx; 
do_plots_check.LABy = InfoLAB.LABy; 
% 'U_ast' = the new reduced basis (U*)
% 'alpha_LS' = coefficients that best reconstruct 'u_mant' (Least-Squares fit)
[U_ast,alpha_LS,velo_LS,norm_dif,contador_svd,final_size] = SVD_ordered_Ubasis(u_mant,reduction,InfoMesh.nel_x,InfoMesh.nel_y,U_basis_norm,InfoMesh,do_plots_check);
disp(['there exist ',num2str(contador_svd),' spare eigenvalues'])
disp(['the relative error in velocities is ',num2str(norm_dif)])
%% PLOT THE ESTIMATION 
% addpath '01_plots'
% varargin.title = 'Temperature estimation';
% plot_results(InfoProblem.T_LAB,1,Temp_est,901,InfoMesh,InfoLAB,varargin);
% daspect([1 1 1])
% title1 = 'Temp. dif. between estimation and reference';
% check_plot_contour(112,Temp_est-Temp,InfoMesh.X,InfoLAB,title1);
% daspect([1 1 1])
%% Calculate error for the Least-Squares (LS) approximation
addpath '04_TempInMesh' 
% Reconstruct the velocity field using the reduced basis (U_ast) and LS coefficients (alpha_LS)
u_LS = reshape(U_ast * alpha_LS, size(u_mant));
% Solve the temperature equation using this approximated velocity
Temp_u_LS = solveTemperatureProblem(u_LS,Temp_est,pres_est,cond_dir,InfoMesh,InfoProblem,InfoMaterial);
% Compute the error (misfit) at the LAB
addpath '06_OptimizationProcedure' 
e_u_LS = compute_errorOnLAB(Temp_u_LS,InfoProblem.T_LAB,InfoMesh.list_cut,0,InfoMesh);
disp(['error function value using velocity approx. with velocity basis: ',num2str(e_u_LS)])
% % plot
% % Least-Squares results (visualize what do they look like)
% varargin.title = 'Results considering input u approximated w/ least-sq';
% plot_results(InfoProblem.T_LAB,1,Temp_u_LS,101,InfoMesh, InfoLAB,varargin);
% daspect([1 1 1])
%% GENERATE INITIAL POINT OF MINIMIZATION
load('98_data_base/variance2comp.mat') % starting guess (alpha_0) for the optimization
ini_data = alpha_LS; % Simple guess: use the LS solution as the start
 
% Calculate the velocity and temp fields for this initial guess
velo_ini_data= reshape(U_ast * ini_data, size(u_mant));
Temp_ini_data = solveTemperatureProblem(velo_ini_data,Temp_est,pres_est,cond_dir,InfoMesh,InfoProblem,InfoMaterial);
%% Plot Initial Guess Results
varargin.title = 'Results considering u from initial \delta data';
varargin.title = '';
addpath '01_plots'
plot_results(InfoProblem.T_LAB,1,Temp_ini_data,108,InfoMesh, InfoLAB,varargin);
daspect([1 1 1])
%% Minimization using Matlab function
addpath '07_Minimization'
nonlinearcase = 1; % 1 = Use non-linear velocity term in optimization
t1 = tic;

% Run the core optimization (e.g., fminunc)
% This finds the optimal 'sol' (alpha coefficients) that minimize the error
% function (misfit at the LAB).
[u_min, sol, fval2, exitflag, output, ini_error, OptName] = run_matlab_minim_vect(nonlinearcase,ini_data, cond_dir, U_ast, InfoMesh, InfoProblem, InfoMaterial,Temp_est,pres_est);
time_mat = toc(t1)
%% Post-Processing: Final Temperature
% Calculate the final temperature field using the optimized velocity 'u_min'
Temp_min_vect = solveTemperatureProblem(u_min,Temp_est,pres_est,cond_dir,InfoMesh,InfoProblem,InfoMaterial);
%% Plot Final Optimization Results
addpath '01_plots'
param.parameters = 0;
param.title = '';
% Plot optimized velocity field
plot_velocities(u_min,1,1012,InfoLAB.LABx,InfoLAB.LABy,InfoMesh.X_v,param);
daspect([1 1 1])
% Plot optimized temperature field
varargin.title = '';
plot_results(InfoProblem.T_LAB,1,Temp_min_vect,1011, InfoMesh,InfoLAB,varargin);
daspect([1 1 1])
% % (Optional) Plot temperature difference
% title1 = 'Temp. dif. between minimization and reference';
%% Calculate and Save Final Errors
diff_u_rel = norm(u_min - u_mant)/norm(u_mant(:));
diff_deltaLS_rel = norm(sol - alpha_LS) /norm(alpha_LS); 
txt1 = ['relative error in alpha Least-squares ',num2str(diff_deltaLS_rel)]; 
txt2 = ['relative error in velocity ',num2str(diff_u_rel)];
disp(txt1)
disp(txt2)