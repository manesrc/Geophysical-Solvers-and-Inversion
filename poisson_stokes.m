function [results,InfoLAB,InfoMesh,InfoProblem] = poisson_stokes(informationLAB,InfoMesh,InfoProblem,T_data,tol1, ...
                                                                                                InfoMaterial,plot_up,InfoMinimization,varargin2)
%[results,InfoLAB,InfoMesh,InfoProblem]

% paths and operations
addpath '00_general_operations'
addpath '01_plots'
addpath '02_isothermInformation'
addpath '03_BulkAndNitscheMatrices'
addpath '04_LAB_operations'
addpath '05_LS_Elem_DOF'
% define "disposition" -> 0: linear, 1: sinusoidal, 2: Data Jeremias plot__
[InfoLAB.LABx,InfoLAB.LABy] = defineLAB(informationLAB);
[InfoLAB.LABx,InfoLAB.LABy] = increase_data_LAB(InfoLAB,50);

%% Definitions for the problem of interest (Section inputs)
L_ref = min([InfoProblem.maxDepth,InfoProblem.cubeSurf]); % [m]

T_sup = T_data.T_sup;                                    % [K]
T_LAB = T_data.T_LAB;                                  % [K]
T_ref = max([T_LAB T_sup 1]);             % [K]
k_1 = T_data.k_1;                                            % [W/(K*m)]     
s_1 = T_data.s_1;                                 % [W/(m^3)]     
s_2 = T_data.s_2;                                 % [W/(m^3)]     

InfoProblem.T_LAB = T_LAB/T_ref;    % LAB definition (fixed temperature)
InfoProblem.T_sup = T_sup/T_ref;    % Fixed superior temperature          
InfoProblem.c = 2;                            % factor for applying Nitsche BC
InfoProblem.T_ref = T_ref;                % reference temperature
InfoProblem.L_ref = L_ref;                % reference 

%% LAB (Section inputs)
dYLAB = abs(1000*mean(InfoLAB.LABy));        % [m]
dTLAB = (InfoProblem.T_LAB-InfoProblem.T_sup)*T_ref;                % [K]
grad_T1 = (dTLAB/dYLAB); % [K/m]
grad_aprox = T_data.grad_aprox; % [K/m]
k_2 = k_1 * grad_T1/grad_aprox;              % flux continuity in LAB [1D idea]
%k_2 = 6;

% reference values to operate adimensionlessly 
k_ref = max(k_1,k_2);
InfoProblem.k_ref = k_ref;
InfoProblem.k1 = k_1/k_ref;  % superior domain conductivity
InfoProblem.k2 = k_2/k_ref; % inferior domain conductivity
% flux
q_flux = k_2*grad_aprox;                                     % [W/m^2]
InfoProblem.q2  = q_flux*((L_ref)/(k_ref*T_ref));   % [dimensionless]
InfoProblem.s1 = s_1*((L_ref^2)/(k_ref*T_ref));     % [dimensionless]
InfoProblem.s2 = s_2*((L_ref^2)/(k_ref*T_ref));     % [dimensionless]

%% The mesh - variables X and T (Section inputs)
% mesh
[InfoMesh.X,InfoMesh.T] = CreaMalla_rectangulo_MTF(InfoMesh);

%% gauss points and shape functions (save structure) (Section inputs)
[MeshIntegration.pospg,MeshIntegration.pespg] = quadrature(InfoMesh.elemType,InfoMesh.nne); 
[MeshIntegration.N,MeshIntegration.Nxi,MeshIntegration.Neta]  = shapeFunctions(InfoMesh.elemType,InfoMesh.nne, MeshIntegration.pospg); 

%% compute the level set of the problem (Section definition of lists)
h_caract = max(abs(InfoMesh.fin_x-InfoMesh.ini_x)/InfoMesh.nel_x,abs(InfoMesh.fin_y-InfoMesh.ini_y)/InfoMesh.nel_y);
InfoProblem.tolerance = h_caract*tol1;
% plot the LS on the mesh:
varargin1.plotLS = 0;  % plot if == 1, not plot else
varargin1.tolPlot = InfoProblem.tolerance;
varargin1.nel_x = InfoMesh.nel_x; 
varargin1.nel_y = InfoMesh.nel_y; 
% Level Set
InfoLAB.maxDepth = L_ref;
InfoMesh.LS_mesh = LevelSet(InfoMesh.X,InfoLAB,varargin1); 

%% find Nitsche elements (Section definition of lists)
plot1 = 0;
[InfoMesh.list1,InfoMesh.list2,InfoMesh.list_cut,InfoMesh.list_edge1,InfoMesh.list_edge2] = findElemNitsche_LS(InfoMesh,InfoLAB,InfoProblem.tolerance,plot1);

%% define dof in sub-domains
[DOF1,DOF2,DOF_interphase,DOF_G] = define_dof(InfoMesh);

%% STEP 1: Solution in Omega_1
problem_int = 1; 
% Obtain FEM matrices
[K_fem1,f_fem1,Ar_phys_1] = K_FEM_n_sparse(InfoMesh,MeshIntegration, InfoProblem,problem_int);
% Nitsche matrix computation
[G1,M1,b1,m1,~] = Nitsche_matrices2D_sparse(InfoMesh, InfoProblem,problem_int,Ar_phys_1);  
N1 = -G1 + M1;
vn1 = -b1+m1;
% FEM + Nitsche
MAT1 = K_fem1+N1;
vect1 = f_fem1 + vn1;
% solve OMEGA_1
t1 = tic;
T_1_nd = MAT1(DOF1,DOF1)\vect1(DOF1);
toc(t1)
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
t2 = tic;
T_2_nd_gmean = MAT2(DOF2,DOF2)\vect2_gmean(DOF2);
toc(t2)
% % dimensional solution OMEGA_1 [in K]
T_2_dim_gmean = T_ref*T_2_nd_gmean; 

%% Stokes problem
plot_up.LAB = [InfoLAB.LABx InfoLAB.LABy];

addpath '06_stokes_2D'
Temp1 = zeros(size(InfoMesh.X,1),1);
Temp2 = zeros(size(InfoMesh.X,1),1);
Temp1(DOF1)  = T_1_dim;
Temp2(DOF2) = T_2_dim_gmean;

plot_up.fig = 0;
plot_up.parameters = 0;

varargin2.u_st = []; varargin2.vel_01= 1;

if (isempty(varargin2.u_st) == 1) && (InfoProblem.vel_01 == 1)
    t_st = tic;
    [u_st,p_st] = ComputeStokesProblem(Temp1,Temp2,InfoMesh,InfoProblem,InfoMaterial,plot_up); % in [m/s] and [Pa]
    toc(t_st)
else
    u_st = varargin2.u_st;
    p_st = varargin2.p_st;
end 

% DISCLAIMER: u_st is originally quadratic, but it's converted to linear before taking it out to use afterwards

%% Advection-diffusion T_2 in Omega_2
u_nd = u_st * ((L_ref*InfoMaterial.rho_ref)/(T_ref*k_ref))^(1/3);
InfoMaterial.calorific = InfoMaterial.calorific_dim * ((L_ref^2 * T_ref * InfoMaterial.rho_ref^2)/ (k_ref^2))^(1/3);

addpath '07_convection_gradT'

if InfoProblem.vel_01 == 1
    % convection matrix:
    G2_matrix = computeConvectionMatrix(u_nd,MeshIntegration,InfoMesh,InfoProblem,InfoMaterial);
    MAT2_new = K_fem2(DOF2,DOF2)+N2(DOF2,DOF2)+G2_matrix(DOF2,DOF2);
else
    MAT2_new = K_fem2(DOF2,DOF2)+N2(DOF2,DOF2);
end 

vect2_indep = f_fem2(DOF2) + vn2(DOF2) + g2_mean(DOF2) ;
T2_2_nd_gmean_adv = MAT2_new\vect2_indep;

%% Compute gradient matrices
% compute gradient matrix such that flux_1 = k_1 * \nabla T_1 = G_1 * T_1
if InfoMesh.nel_x == 100
    num_points_per_elem = 5;
elseif InfoMesh.nel_x == 200
    num_points_per_elem = 2;
elseif InfoMesh.nel_x == 400
    num_points_per_elem = 1;
else 
    num_points_per_elem = 5;
    %warning('DISCLAIMER: this is considered only for doing a convergence (re-consider if not the case)')
end 

addpath '08_Minim_Operations'
figurePoints = 0;       % plot where the points end 
problem_int = 1;    % which matrix to compute
min_ratio = 0.2;     % tolerance of are inside the physical domain note as both domains are 
                              % considered the complementary ratio comp_ratio = 0.8 is also considered
[matG1,PosMat1] = compute_Grad_interphase(num_points_per_elem,DOF1,problem_int,figurePoints,InfoProblem,InfoMesh,Ar_phys_1,min_ratio);
InfoMesh.PosMat = PosMat1;
problem_int = 2;
matG2 = compute_Grad_interphase(num_points_per_elem,DOF2,problem_int,figurePoints,InfoProblem,InfoMesh,Ar_phys_2,min_ratio);
%% just to compare
addpath '10_rest_minim'
results = function_comparisson(MAT2_new,vect2_indep,g2_mean(DOF_G),line_int_vect,matG1,matG2,T_1_nd,T2_2_nd_gmean_adv,DOF1,DOF2,DOF_G,PosMat1,InfoLAB,InfoProblem,InfoMesh,InfoMinimization);

results.g_mean = g2_mean;
results.T1nd = T_1_nd;
results.u_mantle = u_st;
results.p_mantle = p_st;
results.DOF1 = DOF1;
results.DOF2 = DOF2;
results.DOF_interphase = DOF_interphase;
results.DOF_G = DOF_G;
results.PosMat = PosMat1;

