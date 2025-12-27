function G_convMat = computeConvectionMatrix(velo, Temp, pres, InfoMesh,InfoMaterial)
% Creates the matrix computing the convection heat term contracting the velocity with 
% the gradient of the temperature
%
% INPUT:
%   InfoMesh:    X,T
%   InfoMaterial: material characteristics information gauss points + shape functions 
%                  and derivetives
%   InfoProblem: Geometrical characteristics of the problem
%   InfoLAB: Information of the interface
% OUTPUT:
%       G_convMat: Matrix that when multiplied with the temperature applies
%                   the contraction between velocity and gradient of temperature

%% Code
T = InfoMesh.T;
X = InfoMesh.X; 
rho_ref = InfoMaterial.rho_ref;
Tsup = min(Temp);

nnode_el = size(T,2);

% unfold the shape functions
[nelem,nne] = size(T);
nne_v = size(InfoMesh.T_v,2);
elemType = InfoMesh.elemType;
ngp = InfoMesh.ngp_pois;
[pospg,pespg] = quadrature(elemType,ngp);
[N,Nxi,Neta] = shapeFunctions(elemType,nne,pospg);
[N_velo,~,~] = shapeFunctions(elemType,nne_v,pospg);

% % sparse matrices
row_GI = zeros(size(N,2)*size(N,2),nelem);
col_GJ = zeros(size(N,2)*size(N,2),nelem);
val_GX = zeros(size(N,2)*size(N,2),nelem);
% compute convection matrix
for ielem = 1:nelem
    Te = T(ielem,:);
    Te_v = InfoMesh.T_v(ielem,:);
    Xe = X(Te,:);
    velo_element = velo(Te_v,:);
    if (isempty(Temp) && isempty(pres)) == 0
        LS_gp = 0;
        Temp_gp = N * Temp(Te); % [K]
        pres_gp = N*pres(Te);
        rho_gp = DensityAtGaussPoints(LS_gp, Temp_gp,pres_gp,Tsup,InfoMaterial);
    else
        rho_gp = rho_ref * ones(size(N,1),1);
    end 
    Ge = mgElementMatrix(Xe,velo_element, rho_gp, N_velo, N, Nxi, Neta, pespg, nnode_el, InfoMaterial);
    % indices
    [mj,mi] = meshgrid(Te,Te);
    % matrix fill 
    row_GI(:,ielem) = mi(:);
    col_GJ(:,ielem) = mj(:);
    val_GX(:,ielem) = Ge(:); 
end
ndof_t =size(X,1);
G_convMat = sparse(row_GI,col_GJ,val_GX,ndof_t,ndof_t);
end

function Ge = mgElementMatrix(Xe,velocity, rho, N_velocity, N,Nxi,Neta, pespg,numberOfNodes, InfoMaterial)

numberOfGaussPoints = length(pespg); 
Ge = zeros(numberOfNodes, numberOfNodes); 
cal = InfoMaterial.calorific;
    for igaus = 1:numberOfGaussPoints 
        N_igaus = N(igaus,:);
        N_velo_igaus = N_velocity(igaus,:);
        velo_gp = N_velo_igaus * velocity;
        jacob = [Nxi(igaus,:)*Xe(:,1)  Nxi(igaus,:)*Xe(:,2)
            Neta(igaus,:)*Xe(:,1) Neta(igaus,:)*Xe(:,2)]; 
       dvolu = pespg(igaus) * det(jacob); 
       res = jacob\[Nxi(igaus,:); Neta(igaus,:)]; 
       u_gradT = res' * velo_gp';
       % 
       Ge = Ge + rho(igaus) * cal * N_igaus' * u_gradT' *dvolu; 
    end

end
