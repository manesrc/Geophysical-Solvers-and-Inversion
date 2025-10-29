function G_convMat = computeConvectionMatrix(velo,MeshIntegration,InfoMesh,InfoProblem,InfoMaterial)
% Creates the matrix computing the convection heat term contracting the
% velocity with the gradient of the temperature
%
% INPUT:
%   InfoMesh:    X,T
%   MeshIntegration: position and weigth of
% gauss points + shape functions and derivetives
%   InfoProblem: Geometrical characteristics of the problem
%   InfoLAB: Information of the interface
% OUTPUT:
%       G_convMat: Matrix that when multiplied with the temperature applies
%       the contraction between velocity and gradient of temperature

%% Code
T = InfoMesh.T;
X = InfoMesh.X;
LS = InfoMesh.LS_mesh;

nnode_el = size(T,2);


if isempty(InfoMesh.list_cut)~=1 && isempty(InfoMesh.list_edge2) ~= 1
    elementsOmega2 = [InfoMesh.list2; InfoMesh.list_cut(:,1); InfoMesh.list_edge2(:,1)];
    elem_wellCrossed = InfoMesh.list_cut(:,1);
    elem_cross_info = InfoMesh.list_cut;
elseif isempty(InfoMesh.list_cut)
    elementsOmega2 = [InfoMesh.list2; InfoMesh.list_edge2(:,1)];
    elem_wellCrossed = [];
elseif isempty(InfoMesh.list_edge2)
    elementsOmega2 = [InfoMesh.list2; InfoMesh.list_cut(:,1)];
    elem_wellCrossed = InfoMesh.list_cut(:,1);
    elem_cross_info = InfoMesh.list_cut;
else
    error('problems with the elements list')
end

% number of elements
nelem = size(elementsOmega2,1);

% unfold the shape functions
pespg = MeshIntegration.pespg;
N = MeshIntegration.N; 
Nxi = MeshIntegration.Nxi; 
Neta = MeshIntegration.Neta; 

% % sparse matrices
row_GI = zeros(size(N,2)*size(N,2),nelem);
col_GJ = zeros(size(N,2)*size(N,2),nelem);
val_GX = zeros(size(N,2)*size(N,2),nelem);

for ielem1 = 1:nelem
    ielem = elementsOmega2(ielem1);
    Te = T(ielem,:);
    Xe = X(Te,:);
    velo_element = velo(Te,:);
    if ismember(ielem,elem_wellCrossed)        % complicated integration
        interface.elem = 1;  % more GP in the integration
        interface.plot = 0; 
        interface.tryDelaunay = InfoMesh.Delaunay;
        line_elem = ielem == elem_cross_info(:,1);
        interface.P_in_xieta = elem_cross_info(line_elem,2:3); 
        interface.P_out_xieta = elem_cross_info(line_elem,4:5);
        interface.LS_elem = LS(Te); % to assess if the GP are inside the physical part
    else                    % easy integration
        interface.elem = 0;
    end 
    Ge = mgElementMatrix(Xe,velo_element,interface,nnode_el, pespg, N, Nxi, Neta, InfoProblem,InfoMesh,InfoMaterial);
    % indices
   [mj,mi] = meshgrid(Te,Te);
   % matrix fill 
   row_GI(:,ielem1) = mi(:);
   col_GJ(:,ielem1) = mj(:);
   val_GX(:,ielem1) = Ge(:);
end

ndof_t =size(X,1);

G_convMat = sparse(row_GI,col_GJ,val_GX,ndof_t,ndof_t);


end