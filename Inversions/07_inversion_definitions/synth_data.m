%% main inversion code 2
% Define:   
% 1) Domain geometry: lx, ly
% 2) Mesh; 
%%% data:
% 3) Load data to invert
%%% LAB^0:
% 4) initial LAB geometry \GLAB^0
%%% proposal definition
% 5) sigma_d (length of LAB movements); 
% 6) len_prop (length of movement propagation)
%%%% steps
% 7) nsteps

%% FORWARD METHOD DEFINITION
forward_method = 1; % method 1 (split-domain) or method 2 (entire-domain)
%% Code
% domain 
Long_X = 660*1000;
maxDepth  = 660*1000;
% mesh
nel_x = 30;
nel_y = 30;
% load data
%warning('ensure data it is correctly computed')
data = load('data4inversion_linear_nelx30x_nely30.mat');
Data_obs = data.Temp;
% Define initial LAB 

% mx = varargin{1}.mx; 
%         model_height = varargin{1}.model_height;
%         model_width = varargin{1}.model_width; 

depth_ini = maxDepth - 60*1000; 
mx = -0.1;
LAB_definition.rescaleY = 1.1;
LAB_definition.Y_0 = 0; 
LAB_definition.case = 'linear';
LAB_definition.mx = mx;
LAB_definition.depth_ini = depth_ini;

%% inversion definitions
% plots and convection computation
plot_cut_elements = 0; 
velo_01 = 1;
cond_dir = 1; 
% steps
nsteps = 100*1000; 
h = ceil(max(maxDepth/nel_y,Long_X/nel_x));
% proposal definitions
sigma_d = 0.18*h;
len_prop = 2.0*h;
%% Prior definitions
LAB_priorGeometry.LAB_minDepth = 40*1000;        % minimum depth of the LAB
LAB_priorGeometry.LAB_maxDepth = 135*1000;     % maximum depth of the LAB
LAB_priorGeometry.LAB_XChanges = []; % where to input the conditions

%% inversion 1:
inversion_case = 'flux_sup'; % cases: 'Temp_Omega'; 'Temp_sup', 'flux_sup'
make_up = 1;% level set make-up or not (1/0)
[inverted_variable1,error_vect1,phi_t1, accepted_samples1,rejected_samples1,time_inversion1,...
    relevant_inv_dof1, InfoMesh1, InfoProblem1] = inver_function(forward_method,inversion_case, cond_dir, Long_X, maxDepth, nel_x, nel_y, sigma_d, ...
    len_prop, LAB_definition, nsteps, Data_obs, plot_cut_elements, velo_01, LAB_priorGeometry, make_up);
disp(['the ratio of accepted samples is ',num2str(accepted_samples1/nsteps)])

%% save
LAB_definition.LongX = Long_X; 
name2 = generateNames2save(forward_method, inversion_case,len_prop,sigma_d,nel_x,nel_y,nsteps,LAB_definition);
save(name2,'inverted_variable1','error_vect1','phi_t1','accepted_samples1','rejected_samples1','time_inversion1','relevant_inv_dof1','-v7.3')