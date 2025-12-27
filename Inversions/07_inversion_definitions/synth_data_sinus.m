%% FORWARD METHOD DEFINITION
forward_method = 1; % method 1 (split-domain) or method 2 (entire-domain)
%% Code
% domain 
Long_X = 660*1000;
maxDepth  = 660*1000;
% mesh
nel_x = 30;
nel_y = 30;
%LAB
depth_ini = maxDepth - 60*1000; 
mx = -0.1;
LAB_definition.Y_0 = 0; 
LAB_definition.case = 'sinusoidal';
LAB_definition.rescaleY = 1;
LAB_definition.mx = mx;
LAB_definition.model_height = maxDepth;
LAB_definition.model_width = Long_X; 

plot_cut_elements = 0; 
velo_01 = 1;
cond_dir = 1; 

%% load folders and functions
addpath '00_general_operations'
addpath '01_plots'
addpath '02_mesh_cond_SF'
addpath '03_LevelSet'
addpath '04_TempInMesh'
addpath '05_Stokes_2D'
addpath '89_LAB_definition'
addpath '06_Forward1stMethod'
%% inputs
InfoProblem.L_ref = max(Long_X,maxDepth);
% many parameters are in the function below
[InfoProblem, InfoMaterial] = Inputs_file(InfoProblem);
% generate mesh and add more data
InfoMesh = generate_mesh_related_data(nel_x, nel_y,Long_X,maxDepth,InfoProblem);

%% First LAB proposal (m^0)
% compute \Gamma_LAB^0 
InfoLAB0 = defineLAB('sinusoidal',LAB_definition);

grad_T1 = (InfoProblem.T_LAB-InfoProblem.T_sup)/(maxDepth - mean(InfoLAB0.LABy)); % [K/m]
InfoProblem.k2 = InfoProblem.k1 * grad_T1 / InfoProblem.grad_apprT2; % [W/(m K)]        
% compute bc and properties
InfoProblem.q2 = InfoProblem.k2 * InfoProblem.grad_apprT2; % [W/(m^2)]
InfoProblem.T_inf = InfoProblem.T_LAB + InfoProblem.grad_apprT2 * mean(InfoLAB0.LABy);

%% Compute temperature and pressure estimation for matrix computation
varargin.plotLS = 0; varargin.nel_x = InfoMesh.nel_x; varargin.nel_y = InfoMesh.nel_y;
InfoMesh.LS_new = LevelSet(InfoMesh.X,InfoLAB0,varargin);
% find elements crossed by interphase
h_caract = InfoProblem.L_ref * max(abs(InfoMesh.fin_x-InfoMesh.ini_x)/InfoMesh.nel_x,abs(InfoMesh.fin_y-InfoMesh.ini_y)/InfoMesh.nel_y);
InfoProblem.tolerance = h_caract * 0.05;
points2add2LAB = 1000; 
[InfoLAB0.LABx,InfoLAB0.LABy] = increase_data_LAB(InfoLAB0,points2add2LAB);
[InfoMesh.list_Omega1,InfoMesh.list_Omega2,InfoMesh.list_cut, InfoMesh.list_edge1,InfoMesh.list_edge2, interphasePoints] = findElemNitsche_LS_mod1(InfoProblem.tolerance, ...
    plot_cut_elements,InfoMesh,InfoLAB0,InfoProblem);
assert(size(InfoMesh.list_edge1,1) == size(InfoMesh.list_edge2,1))
% estimate temperature and pressure 
[Temp_est,pres_est] = estimate_Temp_pres(InfoMesh.X,InfoMesh.nel_x,InfoMesh.nel_y,InfoMesh.LS_new, InfoMaterial,InfoProblem);



varargin.vel_01 = velo_01;
plot_velo.fig = 0;
plot_velo.LAB = [interphasePoints(1,:); interphasePoints(2,:)] / InfoProblem.L_ref; 
plot_velo.parameters = 0; 

[Temp, results1,InfoMesh1,InfoProblem1] = poisson_stokes_dim(Temp_est, pres_est,InfoMesh,InfoProblem,InfoMaterial,varargin, plot_velo);

%% save
T_Omega1 = results1.T1nd;
T_Omega2 = results1.T2nd_res;
save('data4inversion_sinusoidal_nelx30x_nely30.mat','T_Omega1','T_Omega2','Temp');