function Temp1 = solveTemperatureProblem(u,T_prev,p_prev, cond_dir,InfoMesh,InfoProblem,InfoMaterial)
% Solves the temperature problem in the mesh considering conduction and convection
% OUTPUT:
%   Temp1: updated temperature field [K]
% INPUT:
%   u: velocity field [m/s]
%   T_prev: previous temperature field [K]
%   p_prev: previous pressure field [Pa]
%   cond_dir: direction of the conduction boundary condition (1: bottom and top,
%              2: top dirichlet and bottom neumann ∂T/∂z = 0.5 K/km)
%   InfoMesh: mesh information
%   InfoProblem: problem information
%   InfoMaterial: material information

%% Code
% compute conduction matrix and vector 
[K,f,g_inf] = K_temp_FEM_sparse(T_prev, p_prev, cond_dir, InfoMesh,InfoProblem, InfoMaterial);
fg = f + g_inf;
% compute convection matrix
%assert(size(u,1) == size(InfoMesh.X(:,1),1))


%velo, rho_lin, InfoMesh,InfoMaterial
G_convMatrix = computeConvectionMatrix(u, T_prev, p_prev, InfoMesh,InfoMaterial);
% continue code
all_nodes = 1:1:size(InfoMesh.X,1);
cond_nodes_sup = all_nodes(InfoMesh.X(:,2) == max( InfoMesh.X(:,2) ) );
T_sup = InfoProblem.T_sup; 

if cond_dir == 1
    cond_nodes_inf = all_nodes(InfoMesh.X(:,2) == min(InfoMesh.X(:,2)));
    T_inf = InfoProblem.T_inf; 
    cond_nodes = [cond_nodes_sup cond_nodes_inf];
else
    cond_nodes = cond_nodes_sup;
end 

% Define code to minimize
L_mat = zeros(length(cond_nodes),size(K,2));
cond_value = zeros(length(cond_nodes),1);
for ii=1:length(cond_nodes)
    L_mat(ii,cond_nodes(ii)) = 1;
    if ii > length(cond_nodes_sup)
        cond_value(ii) = T_inf;
    else
        cond_value(ii) = T_sup;
    end 
end 



ZeroMat = zeros(size(L_mat,1));
% Calculate symbolic temperature distribution
Temp1 = [K+G_convMatrix L_mat'; L_mat ZeroMat] \ [fg; cond_value];
Temp1(end-(length(cond_nodes)-1):end) = [];


end 