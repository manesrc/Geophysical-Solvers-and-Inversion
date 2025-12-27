function function_comparation(vel_01,K_fem2,N2,G2_matrix,f_fem2,vn2,matG1,matG2,T_1_nd,DOF1,DOF2,DOF_Q,PosMat,ii_unr,InfoLAB,InfoProblem,InfoMesh)
% create v_0 and B_0
% New Mat 2
MAT2_new2 = K_fem2+N2+vel_01*G2_matrix;
% integrate the flux in the bottom
if InfoMesh.nne == 9
    eleType1D = 0; % line
    nne_1D = 3; % quadratic
elseif InfoMesh.nne == 4
    eleType1D = 0; % line
    nne_1D = 2; % quadratic
end 
[xgp,wgp] = quadrature(eleType1D,nne_1D);
[N,~,~] = shapeFunctions(eleType1D,nne_1D,xgp);
nel_x = InfoMesh.nel_x;
X = InfoMesh.X;
len_x = max(X(:,1)) - min(X(:,1));
x_length_elem = len_x/nel_x; % may be InfoMesh.fin_x-InfoMesh.ini_x instead of 1

% to consider the g_mean integrate over DOF_Q [DOF in model bottom]
line_int_vect = zeros(length(DOF_Q),1);
for jj = 1:nel_x
    if nne_1D == 3
        Te = [2*jj-1 2*jj 2*jj+1];
    else
        Te = [jj jj+1];
    end 
    g_e = zeros(nne_1D,1);
    for kk = 1:size(wgp,2)
        g_e = g_e + N(kk,:)'*wgp(kk) * (x_length_elem/2);
    end
    line_int_vect(Te,1) = line_int_vect(Te,1) + g_e;
end 
g_mean = InfoProblem.q2*line_int_vect;
g_mean_long = [g_mean; zeros(length(DOF2)-length(g_mean),1)];

% vect for v_0 with DOF2 length
vect2_mean = f_fem2(DOF2) + g_mean_long + vn2(DOF2);

fluxOmega1 = matG1 * T_1_nd;
T_2_mean = MAT2_new2(DOF2,DOF2)\vect2_mean;
fluxOmega2_gmean = matG2 * T_2_mean;
% v_0
v_0 = -fluxOmega1 - fluxOmega2_gmean;
% construct B_0
len_q_inf = length(DOF_Q);
D_mat2 = eye(len_q_inf);
short2long = length(DOF2) - len_q_inf;
D_mat = [D_mat2; zeros(short2long,len_q_inf)];
B_0 = matG2 * (MAT2_new2(DOF2,DOF2)\D_mat);
%% Unrestricted LS
% compute the solution unrestricted least-squares:

method = 1;
relationship = 1e-3;
area_cont = [];
[pert_g_unr,~] = obtain_g2(B_0,v_0,area_cont,method,relationship);
pert_g_unr_long = [pert_g_unr; zeros(short2long,1)] ;
vect2_unr = vect2_mean+pert_g_unr_long;
% temperatures in Omega2
T_2_nd_unrest = MAT2_new2(DOF2,DOF2)\vect2_unr;
% fluxes
fluxOmega2_unr = matG2*T_2_nd_unrest;
flux_jump_unr = fluxOmega1 + fluxOmega2_unr;  
%plot: 

gradient_mean = (g_mean./line_int_vect) * ( InfoProblem.T_ref * InfoProblem.k_ref * InfoProblem.k2 / InfoProblem.L_ref ) * (1/(InfoProblem.k_ref * InfoProblem.k2)) * (1000); % [ K / km]
gradient_pert = (pert_g_unr./line_int_vect) * ( InfoProblem.T_ref * InfoProblem.k_ref * InfoProblem.k2 / InfoProblem.L_ref )* (1/(InfoProblem.k_ref * InfoProblem.k2)) * (1000); % [ K / km]


plot_temp_flux(ii_unr,gradient_mean,gradient_pert,T_1_nd,T_2_nd_unrest,flux_jump_unr,fluxOmega1,fluxOmega2_unr,fluxOmega2_gmean,PosMat,DOF1,DOF2,InfoMesh,InfoLAB,InfoProblem);

%% Restricted LS
% compute the solution for restricted least-squares:
Lower_bound = ((0.3-0.5)/0.5)*g_mean;
Upper_bound = ((0.6-0.5)/0.5)*g_mean;
perturbation = optimvar('pert',len_q_inf,'LowerBound',Lower_bound,'UpperBound',Upper_bound);
%residual = (B_0'*B_0+S12(2,2)*D_mat2)*perturbation - B_0'*v_0;
residual = B_0*perturbation-v_0;
obj = residual'*residual;
minim_cond = optimproblem('Objective',obj);
%opts = optimoptions(minim_cond);
%opts.Algorithm = 'trust-region-reflective';
%opts.SubproblemAlgorithm = 'factorization';    
[sol,fval,exitflag,output] = solve(minim_cond);%,'Options',opts);
perturbation_long = D_mat*sol.pert;
vect2_rest = f_fem2(DOF2)  + vn2(DOF2) + g_mean_long + perturbation_long;
% temperatures in Omega2
T_2_nd_rest = MAT2_new2(DOF2,DOF2)\vect2_rest;
% fluxes
fluxOmega2_rest = matG2*T_2_nd_rest;
flux_jump_rest = fluxOmega1 + fluxOmega2_rest;  
%plot: 
ii_r=ii_unr+8;
gradient_pert_analyt = (sol.pert./line_int_vect) * ( InfoProblem.T_ref * InfoProblem.k_ref * InfoProblem.k2 / InfoProblem.L_ref )* (1/(InfoProblem.k_ref * InfoProblem.k2)) * (1000); % [ K / km]
plot_temp_flux(ii_r,gradient_mean,gradient_pert_analyt,T_1_nd,T_2_nd_rest,flux_jump_rest,fluxOmega1,fluxOmega2_rest,fluxOmega2_gmean,PosMat,DOF1,DOF2,InfoMesh,InfoLAB,InfoProblem);

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
    structure1.grad2_minim_analyt = gradient_pert_analyt;
    structure1.flux1 = fluxOmega1;
    structure1.flux2mean = fluxOmega2_gmean;
    structure1.PosMat = PosMat;
    structure1.DOF1 = DOF1;
    structure1.matG2 = matG2;
    structure1.line_int_vect= line_int_vect;
    [likelihood_vect,misfit_vect,acept_perc,due2prior_reject,due2like_reject] = run_MCMC_minimization(nsimu,nDOFQ,prior_g,g0,n_obs, mean_obs, sigma_obs,B_0,v_0,D_mat, MAT2_new2, constant_vector, DOF2,InfoProblem,InfoLAB,InfoMesh,ii_unr, structure1);
end 
