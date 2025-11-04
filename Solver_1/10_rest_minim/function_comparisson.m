function results = function_comparisson(MAT2_new2,vect2_indep,g2_mean,line_int_vect,matG1,matG2,T_1_nd,T_2_gmean,DOF1,DOF2,DOF_G,PosMat,InfoLAB,InfoProblem,InfoMesh,InfoMinimization)
% OUTPUS: a structure with all the relevant results (temperatures, fluxes and bottom gradients)
% INPUTS: MAT2_new2: K_fem + N_nits + G_conv
% vect2_indep: f_fem+m_nits
% g2_mean: the vector g_ref
% line_int_vect: the tributary area of each degree of freedom in Gamma_bot
% matG1: k_1 *gradient * n_1 at the points where the jump of flux is evaluated
% matG2: k_2 *gradient * n_2 at the points where the jump of flux is evaluated
% T_1_nd: solution in Omega_1
% T_2_gmean: solution in Omega2 associated with g2_mean
% DOF1, DOF2, DOF_G: degrees of freedom in Omega1, Omega2 and Gamma_bot
% PosMat: shape functions N evaluated in the points where the jump of flux is evaluated
% InfoLAB, InfoProblem, InfoMesh: general information of the problem
% InfoMinimization: this structure defines which method will be used

% dimensionalizing factors
dimensionalize_flux = (InfoProblem.T_ref * InfoProblem.k_ref) / InfoProblem.L_ref;
dimensionalize_gradient2 = ( dimensionalize_flux/(InfoProblem.k_ref * InfoProblem.k2)) * (1000); % [ K / km], converts g/line_vector to gradient of temperatures in Omega 2

%fluxes 
fluxOmega1_nd = matG1 * T_1_nd;
fluxOmega2_gmean_nd = matG2 * T_2_gmean;
% save fluxes 
results.fluxOmega1_d = fluxOmega1_nd * dimensionalize_flux;
results.fluxOmega2_gmean_d =  fluxOmega2_gmean_nd * dimensionalize_flux;

% B_0 and v_0: 
% v_0
v_0 = -fluxOmega1_nd - fluxOmega2_gmean_nd;
% construct B_0 
len_q_inf = length(DOF_G);
D_matQ = ((InfoMesh.fin_x - InfoMesh.ini_x)/len_q_inf) * eye(len_q_inf);
D_matQ(1,1) = D_matQ(1,1)/2;
D_matQ(len_q_inf,len_q_inf) = D_matQ(len_q_inf,len_q_inf)/2;
%D_matQ = compute_diff_meshes(DOF_G,DOF_Q,InfoMesh);
short2long = length(DOF2) - len_q_inf;
D_mat = [D_matQ; zeros(short2long,size(D_matQ,2))];
B_0 = matG2 * (MAT2_new2\D_mat);

%% Unrestricted LS
% compute the solution unrestricted least-squares:
method = InfoMinimization.method;
relationship = InfoMinimization.relationship;

gradient_mean = (g2_mean./line_int_vect) * dimensionalize_gradient2;
results.gradient_T_mean = gradient_mean;

if method == 4
    [pert_q_rect,pert_q_sq] = obtain_g(B_0,v_0,method,relationship);
    pert_g_rect_nd = D_matQ*pert_q_rect;
    pert_g_sq_nd = D_matQ*pert_q_sq;
    length_extended = length(vect2_indep) - length(pert_g_sq_nd);
    pert_g_unr_long = [pert_g_rect_nd; zeros(length_extended,1)] ;
    vect2_unr = vect2_indep+pert_g_unr_long;
    % rectangular system
    T_2_nd_unrest_rect = MAT2_new2\vect2_unr;
    % fluxes
    fluxOmega2_unr_nd = matG2*T_2_nd_unrest_rect;
    gradient_pert_unr_d = (pert_g_rect_nd./line_int_vect) * dimensionalize_gradient2;

     % square system
    pert_g_unr_long_sq_nd = [pert_g_sq_nd; zeros(length_extended,1)] ;
    vect2_sq_nd = vect2_indep+pert_g_unr_long_sq_nd;
    T_2_nd_unrest_sq = MAT2_new2\vect2_sq_nd;
    % fluxes
    fluxOmega2_sq_nd = matG2*T_2_nd_unrest_sq;     % the flux considers g_mean and minimization
    %flux_jump_unr_sq = fluxOmega1 + fluxOmega2_sq;  
    gradient_pert_sq_d = (pert_g_sq_nd./line_int_vect)  * dimensionalize_gradient2;

    % save results
    results.pert_unr_rect_nd = pert_g_rect_nd;
    results.pert_unr_sq_nd = pert_g_sq_nd;
    % rectangular system
    results.T2nd_unr_rect = T_2_nd_unrest_rect;
    results.gradient_T_pert_unr_rect_d = gradient_pert_unr_d;
    results.fluxOmega2_unr_rect_d = fluxOmega2_unr_nd * dimensionalize_flux; 
    % square system
    results.T2nd_unr_sq = T_2_nd_unrest_sq;
    results.gradient_T_pert_unr_sq_d = gradient_pert_sq_d;
    results.fluxOmega2_unr_sq_d = fluxOmega2_sq_nd * dimensionalize_flux; 
else
    [pert_g_unr_nd,~] = obtain_g(B_0,v_0,method,relationship);
    pert_g_unr_long = [pert_g_unr_nd; zeros(short2long,1)] ;
    vect2_unr = vect2_indep+pert_g_unr_long;
    % temperatures in Omega2
    T_2_nd_unrest = MAT2_new2\vect2_unr;
    % fluxes
    fluxOmega2_unr_nd = matG2*T_2_nd_unrest;
    gradient_pert_unr_d = (pert_g_unr_nd./line_int_vect)  * dimensionalize_gradient2;
    
    % save results
    results.T2nd_unr = T_2_nd_unrest;
    results.pert_unr_nd = pert_g_unr_nd;
    results.gradient_T_pert_unr = gradient_pert_unr_d;
    results.fluxOmega2_unr = fluxOmega2_unr_nd * dimensionalize_flux;
end 



%% Restricted LS
DOF_Q = DOF_G;
% compute the solution for restricted least-squares:
q2_mean = InfoProblem.q2 * ones(size(DOF_G,2),1);       % flux associated with 0.5 K/km
Lower_bound = (0.3-0.5)/(0.5) * q2_mean;       % variation to 0.5 
Upper_bound = (0.6-0.5)/(0.5) * q2_mean;       % flux associated with 0.5 K/km
perturbation = optimvar('pert',length(pert_q_rect),'LowerBound',Lower_bound,'UpperBound',Upper_bound);
residual = B_0*perturbation-v_0;
obj = residual'*residual;
minim_cond = optimproblem('Objective',obj);
%opts = optimoptions(minim_cond);
%opts.Algorithm = 'trust-region-reflective';
%opts.SubproblemAlgorithm = 'factorization';    
[sol,~,~,~] = solve(minim_cond);%,'Options',opts);
perturbation_rest_short = sol.pert;
perturbation_long = D_mat*perturbation_rest_short;
vect2_rest = vect2_indep + perturbation_long;
% temperatures in Omega2
T_2_nd_rest = MAT2_new2\vect2_rest;
% fluxes
fluxOmega2_rest_nd = matG2*T_2_nd_rest;
pert_g_rest_short_nd = D_matQ * sol.pert;
gradient_pert_rest_d = (pert_g_rest_short_nd./line_int_vect) * dimensionalize_gradient2;

% save results
results.T2nd_res = T_2_nd_rest;
results.pert_res = pert_g_rest_short_nd;
results.gradient_T_pert_res_d = gradient_pert_rest_d;
results.fluxOmega2_res_d = fluxOmega2_rest_nd * dimensionalize_flux;

end 