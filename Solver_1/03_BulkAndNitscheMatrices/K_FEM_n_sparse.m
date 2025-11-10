function [K_red,f_red,Ar_phys_n] = K_FEM_n_sparse(InfoMesh, MeshIntegration, InfoProblem,problem_int)
% creates matrix, source term (volumetric: f & Neumann: g) and measures the
% physical part of the element in each subdomain
% INPUT
%   InfoMesh:    X,T
%   MeshIntegration: position and weigth of
% gauss points + shape functions and derivetives
%   InfoProblem: Geometrical characteristics of the problem
%   problem_int: problem in Omega_1 or Omega_2
% OUTPUT
%   K_red:  General matrix with the boundary conditions applied
%   f_red, g_n: indep vector w/ BC and g_n the Neumann BC vector
%   Ar_phys_n: Part of the element physically contained in the domain of interest

%% Code
T = InfoMesh.T;
X = InfoMesh.X;
LS = InfoMesh.LS_mesh;

if problem_int == 1
    Omega_bulk = InfoMesh.list1; % elements completely in Omega_1
    elem_illCrossed = InfoMesh.list_edge1; % elements completely in Omega_1 with interface coinciding with an arista
    elem_wellCrossed = InfoMesh.list_cut; % elements crossed by the interface  
    % Number of elements to compute K_b^1
    nelem = size(Omega_bulk,1) + size(elem_wellCrossed,1) + size(elem_illCrossed,1); 
else
    Omega_bulk = InfoMesh.list2; % elements completely in Omega_2
    elem_illCrossed = InfoMesh.list_edge2; % elements completely in Omega_2 with interface coinciding with an arista
    elem_wellCrossed = InfoMesh.list_cut; % elements crossed by the interface  
    % Number of elements to compute K_b^2
    nelem = size(Omega_bulk,1) + size(elem_wellCrossed,1) + size(elem_illCrossed,1); 
end

if isempty(elem_illCrossed) ~= 1 && isempty(elem_wellCrossed) ~= 1      % elements in both lists
    % elements in the domain
    Omega_int = [Omega_bulk; elem_wellCrossed(:,1); elem_illCrossed(:,1)];
    % store the part of the element (crossed by LAB) contained in the domain 
    Ar_phys_n = zeros(size(elem_wellCrossed,1),2);
    % elements ill-crossed in Omega_1->beta = cte
    list2add = [elem_illCrossed(:,1) ones(size(elem_illCrossed,1),1)];
    Ar_phys_n = [Ar_phys_n; list2add];
elseif isempty(elem_illCrossed) ~= 1 % only ill-crossed elements
    % elements in the domain
    Omega_int = [Omega_bulk; elem_illCrossed(:,1)];
    % elements ill-crossed in Omega_1->Ar_phys = 1 (to be used afterwards)
    list2add = [elem_illCrossed(:,1) ones(size(elem_illCrossed,1),1)];
    Ar_phys_n = list2add;    
    elem_wellCrossed = 0;   % so no element is integrated with more GP
elseif isempty(elem_wellCrossed) ~= 1   % only well-crossed elements
    Omega_int = [Omega_bulk; elem_wellCrossed(:,1)];
    % store the part of the element (crossed by LAB) contained in the domain 
    Ar_phys_n = zeros(size(elem_wellCrossed,1),2);
else
    error('Missing elements')
end 


% unfold the shape functions
pespg = MeshIntegration.pespg;
N = MeshIntegration.N; 
Nxi = MeshIntegration.Nxi; 
Neta = MeshIntegration.Neta; 

% % sparse matrices
row_KI = zeros(size(N,2)*size(N,2),nelem);
col_KJ = zeros(size(N,2)*size(N,2),nelem);
val_KX = zeros(size(N,2)*size(N,2),nelem);

row_fI = zeros(size(N,2),nelem);
col_fJ = ones(size(N,2),nelem);
val_fX = zeros(size(N,2),nelem);

[~,nnode_el] = size(T); 


contador = 0;
for ielem1 = 1:nelem
   ielem = Omega_int(ielem1,1); 
   Te       = T(ielem,:);
   Xe       = X(Te,:);
   if ismember(ielem,elem_wellCrossed(:,1))
       interface.elem = 1;  % more GP in the integration
       interface.plot = 0; 
       interface.tryDelaunay = InfoMesh.Delaunay;
       line_elem = ielem == elem_wellCrossed(:,1);
       interface.P_in_xieta = elem_wellCrossed(line_elem,2:3); 
       interface.P_out_xieta = elem_wellCrossed(line_elem,4:5);
       interface.LS_elem = LS(Te); % to assess if the GP are inside the physical part
       contador = contador+1;
       [Ke,fe,Ar_phys_el]  = mkElementMatrix(Xe, nnode_el, pespg, N, Nxi, Neta, InfoProblem,problem_int,interface,InfoMesh);
       % save the part of the area that is physically in the domain of interest
       Ar_phys_n(contador,:) = [ielem Ar_phys_el];        
   else 
       interface.elem = 0;  % normal integration
       % calculate elemental contribution to stiffness and vectors 
       [Ke,fe,~]  = mkElementMatrix(Xe, nnode_el, pespg, N, Nxi, Neta, InfoProblem,problem_int,interface,InfoMesh);
   end 
   % indices
   [mj,mi] = meshgrid(Te,Te);
   % matrix fill 
   row_KI(:,ielem1) = mi(:);
   col_KJ(:,ielem1) = mj(:);
   val_KX(:,ielem1) = Ke(:);
   % vector fill 
   row_fI(:,ielem1) = Te(:);
   val_fX(:,ielem1) = fe(:);
end 
ndof_t =size(X,1);

K = sparse(row_KI,col_KJ,val_KX,ndof_t,ndof_t);
f = sparse(row_fI,col_fJ,val_fX,ndof_t,1);

all_dof = (1:1:size(K,1));

% boundary conditions in Omega_1
if problem_int == 1
    top_nodes = all_dof(X(:,2) == max(X(:,2))); 
    cond_nodes = ones(size(top_nodes,2),1)*InfoProblem.T_sup;
    DBCmatrix = [top_nodes' cond_nodes];
    [K_red,f_red] = deleteRowsDBC(K,f,DBCmatrix);
else 
    K_red = K;
    f_red = f;
end 

end