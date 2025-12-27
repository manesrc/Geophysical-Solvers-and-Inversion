function results = function_comparation2(MAT2_new2,vect2_mean,g2_mean,line_int_vect,matG1,matG2,T_1_nd,T_2_gmean,DOF1,DOF2,DOF_Q,PosMat,InfoLAB,InfoProblem,InfoMesh,InfoMinimization,plot_resultsTemp)
% create v_0 and B_0
% New Mat 2
%MAT2_new2 = K_fem2+N2+vel_01*G2_matrix;

% vect for v_0 with DOF2 length
%vect2_mean = f_fem2(DOF2) + g2_mean(DOF2) + vn2(DOF2);
fluxOmega1 = matG1 * T_1_nd;
fluxOmega2_gmean = matG2 * T_2_gmean;

% save fluxes 
results.fluxOmega1=fluxOmega1;
results.fluxOmega2_gmean=fluxOmega2_gmean;

% v_0
v_0 = -fluxOmega1 - fluxOmega2_gmean;
% construct B_0
len_q_inf = length(DOF_Q);
D_mat2 = eye(len_q_inf);
short2long = length(DOF2) - len_q_inf;
D_mat = [D_mat2; zeros(short2long,len_q_inf)];

%% Unrestricted LS
% compute the solution unrestricted least-squares:
method = InfoMinimization.method;
relationship = InfoMinimization.relationship;

gradient_mean = (g2_mean./line_int_vect) * ( InfoProblem.T_ref * InfoProblem.k_ref * InfoProblem.k2 / InfoProblem.L_ref ) * (1/(InfoProblem.k_ref * InfoProblem.k2)) * (1000); % [ K / km]
results.gradient_T_mean = gradient_mean;

if method == 4
    [pert_g_rect,pert_g_sq] = obtain_g2(MAT2_new2,D_mat,matG2,v_0,method,relationship);
    pert_g_unr_long = [pert_g_rect; zeros(short2long,1)] ;
    vect2_unr = vect2_mean+pert_g_unr_long;
    % rectangular system
    T_2_nd_unrest_rect = MAT2_new2\vect2_unr;
    % fluxes
    fluxOmega2_unr = matG2*T_2_nd_unrest_rect;
    flux_jump_unr = fluxOmega1 + fluxOmega2_unr;  
    gradient_pert_unr = (pert_g_rect./line_int_vect) * ( InfoProblem.T_ref * InfoProblem.k_ref * InfoProblem.k2 / InfoProblem.L_ref )* (1/(InfoProblem.k_ref * InfoProblem.k2)) * (1000); % [ K / km]

     % square system
     pert_g_unr_long_sq = [pert_g_sq; zeros(short2long,1)] ;
    vect2_sq = vect2_mean+pert_g_unr_long_sq;
    T_2_nd_unrest_sq = MAT2_new2\vect2_sq;
    % fluxes
    fluxOmega2_sq = matG2*T_2_nd_unrest_sq;
    flux_jump_unr_sq = fluxOmega1 + fluxOmega2_sq;  
    gradient_pert_sq = (pert_g_sq./line_int_vect) * ( InfoProblem.T_ref * InfoProblem.k_ref * InfoProblem.k2 / InfoProblem.L_ref )* (1/(InfoProblem.k_ref * InfoProblem.k2)) * (1000); % [ K / km]

    % save results
    results.pert_unr_rect = pert_g_rect;
    results.pert_unr_sq = pert_g_sq;
    % rectangular system
    results.T2nd_unr_rect = T_2_nd_unrest_rect;
    results.gradient_T_pert_unr_rect = gradient_pert_unr;
    results.fluxOmega2_unr_rect = fluxOmega2_unr; 
    % square system
    results.T2nd_unr_sq = T_2_nd_unrest_sq;
    results.gradient_T_pert_unr_sq = gradient_pert_sq;
    results.fluxOmega2_unr_sq = fluxOmega2_sq; 
else
    [pert_g_unr,~] = obtain_g2(MAT2_new2,D_mat,matG2,v_0,method,relationship);
    pert_g_unr_long = [pert_g_unr; zeros(short2long,1)] ;
    vect2_unr = vect2_mean+pert_g_unr_long;
    % temperatures in Omega2
    T_2_nd_unrest = MAT2_new2\vect2_unr;
    % fluxes
    fluxOmega2_unr = matG2*T_2_nd_unrest;
    flux_jump_unr = fluxOmega1 + fluxOmega2_unr;  
    gradient_pert_unr = (pert_g_unr./line_int_vect) * ( InfoProblem.T_ref * InfoProblem.k_ref * InfoProblem.k2 / InfoProblem.L_ref )* (1/(InfoProblem.k_ref * InfoProblem.k2)) * (1000); % [ K / km]
    
    % save results
    results.T2nd_unr = T_2_nd_unrest;
    results.pert_unr = pert_g_unr;
    results.gradient_T_pert_unr = gradient_pert_unr;
    results.fluxOmega2_unr = fluxOmega2_unr;
end 



%% Restricted LS
% compute the solution for restricted least-squares:
Lower_bound = ((0.3-0.5)/0.5)*g2_mean;
Upper_bound = ((0.6-0.5)/0.5)*g2_mean;
perturbation = optimvar('pert',len_q_inf,'LowerBound',Lower_bound,'UpperBound',Upper_bound);
%residual = (B_0'*B_0+S12(2,2)*D_mat2)*perturbation - B_0'*v_0;
B_0 = matG2 * (MAT2_new2\D_mat);
residual = B_0*perturbation-v_0;
obj = residual'*residual;
minim_cond = optimproblem('Objective',obj);
%opts = optimoptions(minim_cond);
%opts.Algorithm = 'trust-region-reflective';
%opts.SubproblemAlgorithm = 'factorization';    
[sol,~,~,~] = solve(minim_cond);%,'Options',opts);
perturbation_rest_short = sol.pert;
perturbation_long = D_mat*perturbation_rest_short;
vect2_rest = vect2_mean + perturbation_long;
% temperatures in Omega2
T_2_nd_rest = MAT2_new2\vect2_rest;
% fluxes
fluxOmega2_rest = matG2*T_2_nd_rest;
flux_jump_rest = fluxOmega1 + fluxOmega2_rest;  
gradient_pert_rest = (perturbation_rest_short./line_int_vect) * ( InfoProblem.T_ref * InfoProblem.k_ref * InfoProblem.k2 / InfoProblem.L_ref )* (1/(InfoProblem.k_ref * InfoProblem.k2)) * (1000); % [ K / km]

% save results
results.T2nd_res = T_2_nd_rest;
results.pert_res = perturbation_rest_short;
results.gradient_T_pert_res = gradient_pert_rest;
results.fluxOmega2_res = fluxOmega2_rest;

%% MCMC
run_MCMC = 0;

if run_MCMC == 1
    constant_vector = f_fem2(DOF2)  + vn2(DOF2) + g_mean_long;
    nsimu = 1e6;
    nDOFQ = length(DOF_Q);
    prior_g = [Lower_bound  Upper_bound];   % where the result of g has to be 
    g0 = (Lower_bound + Upper_bound)/2;
    n_obs = 100;
    mean_obs = 0;
    sigma_obs = 0.05;
    structure1.Plots = 1;
    structure1.Plots_hist = 0 ;
    structure1.T_1_nd = T_1_nd;
    structure1.grad2mean = gradient_mean;
    structure1.grad2_minim_analyt = gradient_pert_rest;
    structure1.flux1 = fluxOmega1;
    structure1.flux2mean = fluxOmega2_gmean;
    structure1.PosMat = PosMat;
    structure1.DOF1 = DOF1;
    structure1.matG2 = matG2;
    structure1.line_int_vect= line_int_vect;
    [likelihood_vect,misfit_vect,acept_perc,due2prior_reject,due2like_reject] = run_MCMC_minimization(nsimu,nDOFQ,prior_g,g0,n_obs, mean_obs, sigma_obs,B_0,v_0,D_mat, MAT2_new2, constant_vector, DOF2,InfoProblem,InfoLAB,InfoMesh,ii_unr, structure1);
end 

%% plot results? 

if plot_resultsTemp == 1
    %plot: 
    plot_temp_flux(ii_unr,gradient_mean,gradient_pert,T_1_nd,T_2_nd_unrest,flux_jump_unr,fluxOmega1,fluxOmega2_unr,fluxOmega2_gmean,PosMat,DOF1,DOF2,InfoMesh,InfoLAB,InfoProblem);
    %plot: 
    plot_temp_flux2(ii_r,gradient_mean,gradient_pert_rest,T_1_nd,T_2_nd_rest,flux_jump_rest,fluxOmega1,fluxOmega2_rest,fluxOmega2_gmean,PosMat,DOF1,DOF2,InfoMesh,InfoLAB,InfoProblem);
end 


end 