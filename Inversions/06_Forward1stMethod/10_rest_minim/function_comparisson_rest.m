function results = function_comparisson_rest(MAT2_new2,vect2_indep,g2_mean,line_int_vect,matG1,matG2,T_1_nd,T_2_gmean,DOF2,DOF_G,InfoProblem,InfoMesh)
% OUTPUS: a structure with all the relevant results (temperatures, fluxes and bottom gradients)
% INPUTS: MAT2_new2: K_fem + N_nits + G_conv
% vect2_indep: f_fem+m_nits
% g2_mean: the vector g_ref
% line_int_vect: the tributary area of each degree of freedom in Gamma_bot
% matG1: k_1 *gradient * n_1 at the points where the jump of flux is evaluated
% matG2: k_2 *gradient * n_2 at the points where the jump of flux is evaluated
% T_1_nd: solution in Omega_1
% T_2_gmean: solution considering g2_mean
% DOF1, DOF2, DOF_G: degrees of freedom in Omega1, Omega2 and Gamma_bot
% PosMat: N evaluated in the points where the jump of flux is evaluated
% InfoLAB, InfoProblem, InfoMesh: general information of the problem
% InfoMinimization: this structure defines which method will be used


%fluxes 
fluxOmega1_nd = matG1 * T_1_nd;
fluxOmega2_gmean_nd = matG2 * T_2_gmean;
% save fluxes 
dimensionalize_flux = 1;
results.fluxOmega1_d = fluxOmega1_nd * dimensionalize_flux;
results.fluxOmega2_gmean_d =  fluxOmega2_gmean_nd * dimensionalize_flux;

% B_0 and v_0: 
% v_0
v_0 = -fluxOmega1_nd - fluxOmega2_gmean_nd;
% construct B_0
len_q_inf = length(DOF_G);
DOF_Q = 1:1:(InfoMesh.nel_x+1);
D_matQ = (InfoProblem.L_ref) * compute_diff_meshes(DOF_G,DOF_Q,InfoMesh);
short2long = length(DOF2) - len_q_inf;
D_mat = [D_matQ; zeros(short2long,size(D_matQ,2))];
B_0 = matG2 * (MAT2_new2\D_mat);

k2 = InfoProblem.k2;    % [W / (m*K)]
dimensionalize_gradient2 = 1;
gradient_mean = (g2_mean./(k2 * line_int_vect)) * dimensionalize_gradient2; % [K/m]
results.gradient_T_mean = gradient_mean;

%% Restricted LS
% compute the solution for restricted least-squares:
q2_mean = (InfoProblem.q2 * ones(size(DOF_Q,2),1));       % flux associated with 0.5 K/km
Lower_bound = ((0.3-0.5)/0.5) * q2_mean;       % variation to 0.5 
Upper_bound = ((0.6-0.5)/0.5)* q2_mean;       % flux associated with 0.5 K/km
perturbation = optimvar('pert',size(B_0,2),'LowerBound',Lower_bound,'UpperBound',Upper_bound);
residual = B_0*perturbation-v_0;
obj = residual'*residual;
minim_cond = optimproblem('Objective',obj);
% opts = optimoptions(minim_cond);
% opts.Algorithm = 'trust-region-reflective';
% opts.SubproblemAlgorithm = 'factorization';    
%[sol,~,~,~] = solve(minim_cond);%,'Options',opts);

% Define solver options to suppress output
options = optimoptions('quadprog', 'Display', 'none'); % Assuming 'quadprog' is the underlying solver

% Solve the optimization problem with the specified options
warning off
[sol,~,~,~] = solve(minim_cond, 'Options', options);
warning on
perturbation_rest_short = sol.pert;
perturbation_long = D_mat*perturbation_rest_short;
vect2_rest = vect2_indep + perturbation_long;
% temperatures in Omega2
T_2_nd_rest = MAT2_new2\vect2_rest;
% fluxes
fluxOmega2_rest_nd = matG2*T_2_nd_rest;
pert_g_rest_short_nd = D_matQ * sol.pert;
gradient_pert_rest_d = (pert_g_rest_short_nd./(k2*line_int_vect)) * dimensionalize_gradient2;

% save results
results.T2nd_res = T_2_nd_rest;
results.pert_res = pert_g_rest_short_nd;
results.gradient_T_pert_res_d = gradient_pert_rest_d;
results.fluxOmega2_res_d = fluxOmega2_rest_nd * dimensionalize_flux;

end 