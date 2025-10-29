% generate a function that runs the given case and saves the results
%% define important inputs for thermal problem
% unchanging parameters in most of the runs: 
T_data.T_sup = 293; 
T_data.T_LAB = 1573;
T_data.k_1 = 1; 
T_data.s_1 = 0; 
T_data.s_2 = 0;
T_data.grad_aprox = 0.5/1000;   % [K/m]
InfoProblem.maxDepth = 660*1000;
InfoProblem.cubeSurf = 660*1000;
InfoMesh.elemType = 1; InfoMesh.nne = 4; % Linear quads
InfoMesh.ini_x = 0; InfoMesh.fin_x = 1; InfoMesh.ini_y = 0; InfoMesh.fin_y = 1;
InfoMesh.ngp_all = 9; % number of gauss points per element
% interface elements:
InfoMesh.Delaunay = 1;
InfoMesh.ngp_Nits = 9; % number of gauss points in elements crossed by LAB [Bulk terms]
InfoMesh.elemType_LAB = 0; % type of element in the interface LAB [dimension: d-1, means line for 2D]
InfoMesh.ngp_interface = 9; % number of gauss points in interface LAB [dimension: d-1, means line for 2D]
% Stokes problem
plot_up.fig = 0; % plots the quiver with mantle velocities
plot_up.parameters = 0; % plots the paramater value distribution 
InfoMaterial.rho_ref = 3300; % [kg/m^3]
InfoMaterial.mu_ref = 1e24; % [Pa*s]
InfoMaterial.gravity_units = -9.81; % [m/s^2]
InfoMaterial.mu1 = 1e25;
InfoMaterial.mu2 = 1e19;
InfoMaterial.press_ref = 0; % [Pa] = [N/m^2]
InfoMaterial.alpha = 1e-5;
InfoMaterial.beta = 1e-5;
InfoMaterial.mu_max = 1e25;
InfoMaterial.disloc = 1.1e5;
InfoMaterial.second_inv = 1e-15;
InfoMaterial.energy = 5.3e5;
InfoMaterial.volume = 14;   
InfoMaterial.gas_constant = 8.314;
InfoMaterial.flow_law = 3.5;
InfoMaterial.T_ref0 = 293; % [K];
InfoMaterial.calorific_dim = 1200; % [J/(kg*K)]

InfoProblem.vel_01 = 1; % consider convection or not

% to be changed or note
tol1 = 0.05; % considers if the element is too badly cut and disregards it

% plot yes/no
plot_up.temp_plots = 0; % to see plots use 1, case 0 is better for server runs

% LAB definition
informationLAB.disposition = 1; % disposition = 0 (linear), define first (y_ini) and last point (y_fin) [ (1) Sinusoidal; (2) Andean-Jeremias]
informationLAB.y_ini = 0.1; informationLAB.y_fin = 0.2; 
% Mesh and minimization 
InfoMesh.NewMesh = []; % this is to do a transformation from linear to quadratic mesh to pass LBB condition when solving Stokes
InfoMinimization.method = 4; %     cases: method =va 1 rectangular w/ K inverse B * p = v , eigval permuted to "InfoMinimization.relationship"
                                               %                method = 2 square system with (B^T * B + alpha * I)* p = B^Tv,  alpha = relationship^2 * eig_max                                                
                                               %                NOT AVAILABLE method = 3 rectangular w/ G_2 pseudo inverse [NO VA BIEN]
                                               %                method = 4 solves both 1 and 2
InfoMinimization.relationship = 1e-2;       % tolerance for eigenvalues to be considered 

%% SOLVE case for different meshes
% % Mesh 100:
InfoMesh.nel_x = 100; % number of elements in X direction
InfoMesh.nel_y = 100; % number of elements in Y direction

%%
addpath 00_general_operations

if exist('results_mesh_100.mat','file') ~= 0 
     old_result = loadStructFromFile('results_mesh_100','results1');
     varargin.u_st = old_result.u_mantle;
     varargin.p_st = old_result.p_mantle;
elseif exist('results3.mat','file') == 1
    %t0 = tic();
    [u_mant1,p_mant1] = FineMantleVelo2Coarse(results3,InfoMesh3,InfoMesh1);
    %toc(t0)
    varargin.u_st = u_mant1;
    varargin.p_st = p_mant1;
else
    varargin.u_st = [];
    varargin.p_st = [];
  
end 



%%

% run code
[results1,InfoLAB1,InfoMesh1,InfoProblem1] = poisson_stokes(informationLAB,InfoMesh,InfoProblem,T_data,tol1,InfoMaterial,plot_up,InfoMinimization,varargin);
% save results
if InfoProblem.vel_01 == 1
    save 'results_mesh_100' InfoLAB1 InfoMesh1 InfoProblem1 results1
else
    save 'results_mesh_100_no_conv' InfoLAB1 InfoMesh1 InfoProblem1 results1
end 


%%  Post -  Process plots
 
%plot mesh 100:
%RECTANGULAR
perturbation_rect1 = results1.gradient_T_pert_unr_rect_d;
T2_rect1 = results1.T2nd_unr_rect;
fluxOmega2_rect1 = results1.fluxOmega2_unr_rect_d;
ii_figure_num = 10; % this number is the one used to print the results in figures çvarargin.restricted = 0; 
varargin.restricted = 0; % only to fix the plot scale 
plot_temp_flux(ii_figure_num,T2_rect1,perturbation_rect1,fluxOmega2_rect1,results1,InfoMesh1,InfoLAB1,InfoProblem1,varargin); 
%SQUARE 
perturbation_sq1 = results1.gradient_T_pert_unr_sq_d;
T2_sq1 = results1.T2nd_unr_sq;
fluxOmega2_sq1 = results1.fluxOmega2_unr_sq_d; 
varargin.restricted = 0; 
ii_figure_num = ii_figure_num + 50;
plot_temp_flux(ii_figure_num,T2_sq1,perturbation_sq1,fluxOmega2_sq1,results1,InfoMesh1,InfoLAB1,InfoProblem1,varargin);
%RESTRICTED
perturbation_res1 = results1.gradient_T_pert_res_d;
T2_res1 = results1.T2nd_res;
fluxOmega2_res1 = results1.fluxOmega2_res_d;
ii_figure_num = ii_figure_num + 50;
varargin.restricted = 1; 
plot_temp_flux(ii_figure_num,T2_res1,perturbation_res1,fluxOmega2_res1,results1,InfoMesh1,InfoLAB1,InfoProblem1,varargin);
