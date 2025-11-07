function [K,G,f,mu_gp_all,rho_gp_all,xpg_all] = makeGlobalMatrices(T_prev,p_prev, ...
    elemType,X,T,XP,TP,LS_mesh_unitless,nGPoiss,plot_up, material)
% makeGlobalMatrices: assembles global matrices for the 2D Stokes problem
% Outputs:
%   K: global velocity-velocity matrix
%   G: global velocity-pressure matrix
%   f: global force vector
%   mu_gp_all: dynamic viscosity at all gauss points (for plotting)
%   rho_gp_all: density at all gauss points (for plotting)
%   xpg_all: coordinates of all gauss points (for plotting)
% Inputs:
%   T_prev: previous temperature field
%   p_prev: previous pressure field
%   elemType: type of finite element
%   X: coordinates of the mesh nodes
%   T: connectivity matrix for velocity nodes
%   XP: coordinates of the pressure nodes
%   TP: connectivity matrix for pressure nodes
%   LS_mesh_unitless: level set function values at pressure nodes
%   nGPoiss: number of Gauss points for Poisson integration
%   plot_up: structure with plotting options
%   material: structure with material properties

% number of space dimensions
nsd = size(X,2);

% number of nodes
nOfNodes = size(X,1);
[nOfElements,nOfVelocityNodesPerElement] = size(T);
[~,nOfPressureNodesPerElement] = size(TP);

% number of degrees of freedom
nOfVelocityDofPerElement = nsd*nOfVelocityNodesPerElement;

% element matrix size
mKe = nOfVelocityDofPerElement;
nKe = nOfVelocityDofPerElement;

mGe = nOfPressureNodesPerElement;
nGe = nOfVelocityDofPerElement;

% initialisation
allK = zeros(mKe*nKe, nOfElements);
allKI = zeros(mKe*nKe, nOfElements);
allKJ = zeros(mKe*nKe, nOfElements);

allG = zeros(mGe*nGe, nOfElements);
allGI = zeros(mGe*nGe, nOfElements);
allGJ = zeros(mGe*nGe, nOfElements);

f = zeros(nOfNodes*nsd,1);

mu_gp_all = [];
rho_gp_all = [];
xpg_all = [];

% to define: discontinuous depending on LS or integrate nodal values with shape functions
mu_ref = material.mu_ref;
rho_ref = material.rho_ref;
T_sup = min(T_prev);

% loop in standard elements
for ielem = 1:nOfElements
   % velocity nodes 
   Te = reshape([2*T(ielem,:)-1; 2*T(ielem,:)],1,nOfVelocityDofPerElement);
   % pressure nodes 
   TeP = TP(ielem,:);
   % coords of the velocity nodes
   Xe = X(T(ielem,:),:);
   LS_elem = LS_mesh_unitless(TeP);
      
   % number of materials in the element
    nOfMaterialsInElement = [abs(sum(sign(LS_elem(LS_elem~=0)))) == sum(LS_elem ~= 0) ...
        abs(sum(sign(LS_elem(LS_elem~=0)))) ~= sum(LS_elem ~= 0)] * [1; 2]; 
    
    % gauss points and shape functions
    if nOfMaterialsInElement == 1
        % Number of integration points
        nGaussPoints = nGPoiss;
        % quadrature and shape functions
        [pospg1,pespg1] = quadrature(elemType,nGaussPoints);
        [N1,Nxi1,Neta1] = shapeFunctions(elemType,nOfVelocityNodesPerElement,pospg1);
        [NP1,~,~] = shapeFunctions(elemType,nOfPressureNodesPerElement,pospg1);
    else
        % Add Gauss points for a more detailed description
        nGaussPoints = 25;
        % quadrature and shape functions
        [pospg1,pespg1] = quadrature(elemType,nGaussPoints);
        [N1,Nxi1,Neta1] = shapeFunctions(elemType,nOfVelocityNodesPerElement,pospg1);
        [NP1,~,~] = shapeFunctions(elemType,nOfPressureNodesPerElement,pospg1);   
    end 
    % material dependence of temperature and pressure check 
    if (material.rho_pres_temp == 1) || (material.mu_pres_temp == 1)
        % Element and pressure 
        T_elem = T_prev(TeP);
        p_elem = p_prev(TeP);
        Temp_gp = NP1 * T_elem;
        pres_gp = NP1 * p_elem;
    else
        Temp_gp = [];
        pres_gp = []; 
    end 
    LS_gp = NP1 * LS_elem;
    if isfield(material,'velo_basis')
        material.Xgp = NP1*Xe(1:4,:);
    end 
    mu_gp = ViscosityAtGaussPoint(LS_gp,Temp_gp,pres_gp,material);
    rho_gp = DensityAtGaussPoints(LS_gp, Temp_gp,pres_gp, T_sup, material);
    mu_adim = mu_gp / mu_ref; rho_adim = rho_gp / rho_ref;
   % save parameter values for plottingN
  if plot_up.parameters == 1
      mu_gp_all = [mu_gp_all; mu_adim];
      rho_gp_all = [rho_gp_all; rho_adim];
      xpg_all = [xpg_all; N1*Xe];
   else
      mu_gp_all = [];
      rho_gp_all = [];
      xpg_all = [];
  end 
  
   % element matrix
   gravity = material.gravity; 
   [Ke,Ge,fe] = makeElementMatrix(Xe, ...
      nOfVelocityDofPerElement,pospg1,pespg1,N1,Nxi1,Neta1,...
      nOfPressureNodesPerElement,NP1,gravity,rho_adim,mu_adim);

   % store elemental matrices 
   allK(:,ielem) = Ke(:);
   [mi,mj] = meshgrid(Te, Te);
   allKI(:,ielem) = mi(:);
   allKJ(:,ielem) = mj(:);
    
   allG(:,ielem) = Ge(:);
   [mi,mj] = meshgrid(Te, TeP);
   allGI(:,ielem) = mi(:);
   allGJ(:,ielem) = mj(:);

   % assembly force term
   f(Te) = f(Te) + fe;
end

% assembly matrices
K = sparse(allKJ, allKI, allK);
G = sparse(allGJ, allGI, allG);

