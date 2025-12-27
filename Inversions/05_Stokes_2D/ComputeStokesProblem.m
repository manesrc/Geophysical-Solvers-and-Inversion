function [velo_quad,pres] = ComputeStokesProblem(T_est, p_est, InfoMesh,InfoProblem,InfoMaterial,plot_up)
% ComputeStokesProblem: solves the Stokes problem for a 2D domain to compute 
% velocity of the mantle flow and the pressure field.
% OUTPUTS: 
% u: velocity as a vector of DOF in mesh x 2 directions [nDOF_X,2]
% p: pressure as a vector of DOF in mesh x 1 [nDOF_X,1]
% INPUTS:
% LS: The level set values to decide whether the material properties are 1 or 2
% InfoMesh: characteristics of the mesh, the mesh, etc
% InfoProblem: characteristics of the problem, lengths 
% InfoMaterial: description to obtain material characteristics such to compute 
%               rho and mu, etc

%% Meshes
% mesh X has units therefore to solve Stokes needs to be unitless
L_ref = InfoProblem.L_ref;
XP = InfoMesh.X / L_ref;
TP = InfoMesh.T;
X = InfoMesh.X_v / L_ref;
T = InfoMesh.T_v;

%% material properties
%mu_ref = max(mu_lin); % viscosity 
%rho_ref = min(rho_lin); % density [top_domain; bottom_domain] 
mu_ref = InfoMaterial.mu_ref;
rho_ref = InfoMaterial.rho_ref;
InfoMaterial.gravity = InfoMaterial.gravity_units*( (L_ref^3) * (rho_ref^2)/ (mu_ref^2) ); % gravity

%% Build matrices
elemType = InfoMesh.elemType; 
nGPoiss = InfoMesh.ngp_pois;
LS_mesh_temp = InfoMesh.LS_new;
[K,G,f,plot_up.mu_used,plot_up.rho_used,plot_up.xgp22] = makeGlobalMatrices(T_est, p_est, ...
    elemType, X,T,XP,TP,LS_mesh_temp,nGPoiss,plot_up,InfoMaterial);

%% Boundary conditions
ksize = size(K,1);
[Adir,bdir] = bcFreeSlip(X,ksize,XP);
% build the global system
gsize = size(G,1);
Z01 = sparse(gsize,gsize);
Z0cc = sparse( size(Adir,1),size(Adir,1) );

Atot = [[K G'; G Z01; Adir] [Adir'; Z0cc] ];
btot = [f; zeros(gsize,1); bdir];

% solve it!
aux = Atot\btot;

% extract the velocity 
velo1 = aux(1:ksize) * (mu_ref/(L_ref*rho_ref)); % -> back to [m/s], [Pa*s]/([m] * [kg/m^3])
velo_quad = reshape(velo1,2,[])';

% extract pressures
pres = aux(ksize+1:ksize+gsize) * ( (mu_ref^2) / ( (L_ref^2)*rho_ref )  );  % -> back to [Pa]
pres = pres*(1e-6); % [MPa]

%% Post process
if plot_up.fig_velo == 1
    PropChangeX = plot_up.LAB(:,1)*InfoProblem.L_ref;
    PropChangeY = plot_up.LAB(:,2)*InfoProblem.L_ref;
    plot_up.rho_ref = InfoMaterial.rho_ref;
    plot_up.mu_ref = InfoMaterial.mu_ref;
    num_fig = 89;
    addpath '01_plots'
    sc = 1;
    plot_velocities(velo_quad,sc,num_fig,PropChangeX,PropChangeY,X*InfoProblem.L_ref,plot_up)
    figure(num_fig)
end



end