function [K,f,g_vect] = K_temp_FEM_sparse(Temp, pres, cond_dir,InfoMesh,InfoProblem, InfoMaterial)
% creates matrix, source term (volumetric: f & Neumann: g) and measures the
% physical part of the element in each subdomain
% 
% INPUT
%   InfoMesh:    X,T
%   InfoProblem: Geometrical characteristics of the problem
%   InfoLAB: Information of the interface
%
% OUTPUT
%   K_red:  General matrix with the boundary conditions applied
%   f_red, g_n: indep vector w/ BC and g_n the Neumann BC vector
%   Ar_phys_n: Part of the element physically contained in the domain of interest

%% Code
T = InfoMesh.T;
X = InfoMesh.X;
LS = InfoMesh.LS_new;

% unfold the shape functions
[nelem, nne] = size(T);
elemType = InfoMesh.elemType;
ngp = InfoMesh.ngp_pois;
[pospg,pespg] = quadrature(elemType,ngp);
[N,Nxi,Neta] = shapeFunctions(elemType,nne,pospg);
% sparse matrices
row_KI = zeros(size(N,2)*size(N,2),nelem);
col_KJ = zeros(size(N,2)*size(N,2),nelem);
val_KX = zeros(size(N,2)*size(N,2),nelem);
row_fI = zeros(size(N,2),nelem);
col_fJ = ones(size(N,2),nelem);
val_fX = zeros(size(N,2),nelem);

% cut elements
interface = InfoMesh.list_cut; 

k1 = InfoProblem.k1;
k2 = InfoProblem.k2; 

for ielem = 1:nelem
   Te       = T(ielem,:);
   Xe       = X(Te,:);
   LS_elem = LS(Te);
   one_material = (length(Te) == abs(sum(sign(LS_elem))) );
   if ismember(ielem,interface(:,1))
       one_material = 0; 
       %assert();
       relev_data = interface(ielem == interface(:,1),:);   % contains: [elem_num, xi_entry, eta_entry, xi_out, eta_out, n_x, ny]
   else
       one_material = 1; 
       relev_data = [];
   end 
   % temperature and pressure 
   LS_gp = N * LS_elem;
   conductivity = compute_k_gp(one_material,relev_data,k1,k2,LS_gp,Temp, pres,N,Te,InfoMesh.nel_x,InfoMaterial);
   assert(InfoProblem.s1 == InfoProblem.s2)
   source = InfoProblem.s1;
    % calculate elemental contribution to stiffness and vectors 
   [Ke,fe]  = mkElementMatrix(Xe, conductivity, source, nne, pespg, N, Nxi, Neta);
   % indices
   [mj,mi] = meshgrid(Te,Te);
   % matrix fill 
   row_KI(:,ielem) = mi(:);
   col_KJ(:,ielem) = mj(:);
   val_KX(:,ielem) = Ke(:);
   % vector fill 
   row_fI(:,ielem) = Te(:);
   val_fX(:,ielem) = fe(:);
end 
ndof_t =size(X,1);
K = sparse(row_KI,col_KJ,val_KX,ndof_t,ndof_t);
f = sparse(row_fI,col_fJ,val_fX,ndof_t,1);
all_dof = (1:1:size(K,1));
DOF_G = all_dof( X(:,2)==min(X(:,2)) );
if cond_dir == 0 
    g_vect2 = compute_g_mean(DOF_G',InfoMesh,InfoProblem);
    g_vect = sparse(size(f,1),size(f,2));
    g_vect(DOF_G) = g_vect2;
else
    g_vect = sparse(size(f,1),size(f,2));
end 


end

function [Ke,fe] = mkElementMatrix(Xe, conductivity, source, numberOfNodes,pespg, N, Nxi, Neta)
%% elemental matrix calculation
%Xe, nnode, pospg, pespg, N, Nxi, Neta, InfoProblem,problem_int,interface,InfoLAB
numberOfGaussPoints = length(pespg); 
Ke = zeros(numberOfNodes, numberOfNodes); 
fe = zeros(numberOfNodes, 1);
for igaus = 1:numberOfGaussPoints 
    jacob = [Nxi(igaus,:)*Xe(:,1)  Nxi(igaus,:)*Xe(:,2)
            Neta(igaus,:)*Xe(:,1) Neta(igaus,:)*Xe(:,2)]; 
    dvolu = pespg(igaus) * det(jacob); 
    res = jacob\[Nxi(igaus,:); Neta(igaus,:)]; 
    Nx = res(1,:); 
    Ny = res(2,:);
    %
    Ke = Ke + conductivity(igaus)*(Nx'*Nx + Ny'*Ny)*dvolu; 
    fe = fe + source * N(igaus,:)' *dvolu; 
end 
end