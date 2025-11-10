function [K,G,f,mu,temp,rho,xgp22] = makeGlobalMatrices(elemType,X,T,XP,TP,material,T_Omega1,T_Omega2,LS_mesh,InfoMesh,InfoProblem,plot_up)
% makeGlobalMatrices: assembles global matrices for the 2D Stokes problem
%
% INPUT:
% elemType                  type of element
% X                         coordinates of nodes
% T                         connectivity matrix for velocity nodes
% XP                        connectivity matrix for pressure nodes
% TP                        connectivity matrix for pressure nodes
% material                  structure containing material properties
% T_Omega1                  nodal temperature values in subdomain Omega1
% T_Omega2                  nodal temperature values in subdomain Omega2
% LS_mesh                   level set values at mesh nodes
% InfoMesh                  structure containing mesh information
% InfoProblem               structure containing problem information
% plot_up                   flag for plotting updates
%
% OUTPUT:
% K                         global stiffness matrix
% G                         global gradient matrix 
% f                         global force vector
% mu                        viscosity values at gauss points (for postprocessing)
% temp                      temperature values at gauss points (for postprocessing)
% rho                       density values at gauss points (for postprocessing)
% xgp22                     coordinates of gauss points (for postprocessing)    

% number of space dimensions
nsd = 2;

% number of nodes
nOfNodes = size(X,1);
[nOfElements,nOfVelocityNodesPerElement] = size(T);
[~,nOfPressureNodesPerElement] = size(TP);

% number of degrees of freedom

nOfVelocityDofPerElement = nsd*nOfVelocityNodesPerElement;

% Number of integration points
nGaussPoints = 9;

% quadrature and shape functions
[pospg,pespg] = quadrature(elemType,nGaussPoints);
[N,Nxi,Neta] = shapeFunctions(elemType,nOfVelocityNodesPerElement,pospg);
[NP,~,~] = shapeFunctions(elemType,nOfPressureNodesPerElement,pospg);

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

% elements in Omega1 
if (isempty(InfoMesh.list_edge1) ~=1) && (isempty(InfoMesh.list_cut) ~= 1)  % both well and badly cut by the interface 
    listOmega1 = [InfoMesh.list1; InfoMesh.list_edge1(:,1)];
    listOmega2 = [InfoMesh.list2; InfoMesh.list_edge2(:,1)];
    listInt = InfoMesh.list_cut(:,1);
elseif isempty(InfoMesh.list_cut) ~= 1      % all elements are well cut by the interface 
    listOmega1 = InfoMesh.list1;
    listOmega2 = InfoMesh.list2;
    listInt = InfoMesh.list_cut(:,1);
elseif isempty(InfoMesh.list_edge1) ~= 1        % only badly cut
    listOmega1 = [InfoMesh.list1; InfoMesh.list_edge1(:,1)];
    listOmega2 = [InfoMesh.list2; InfoMesh.list_edge2(:,1)];
    listInt = [];
else
    error('missing something')
end 

mu = [];
rho = [];
temp = [];
xgp22 = [];

% loop in standard elements
for ielem = 1:nOfElements
   % velocity nodes 
   Te = reshape([2*T(ielem,:)-1; 2*T(ielem,:)],1,nOfVelocityDofPerElement);
   % pressure nodes 
   TeP = TP(ielem,:);
   % coords of the velocity nodes
   Xe = X(T(ielem,:),:);
   % temperature
   if ismember(ielem,listOmega1)
       temperature_el = T_Omega1(T(ielem,:));
       LS_el = 1;
       tolerance = InfoProblem.tolerance;

                 % element matrix
       [Ke,Ge,fe,muee,tempee,rhoee,Xee] = makeElementMatrix(Xe, ...
          nOfVelocityDofPerElement,pospg,pespg,N,Nxi,Neta,...
          nOfPressureNodesPerElement,NP,material,temperature_el,LS_el,tolerance);
   elseif ismember(ielem,listOmega2)
       temperature_el = T_Omega2(T(ielem,:));
       LS_el = 2;
       tolerance = InfoProblem.tolerance;
                 % element matrix
       [Ke,Ge,fe,muee,tempee,rhoee,Xee] = makeElementMatrix(Xe, ...
          nOfVelocityDofPerElement,pospg,pespg,N,Nxi,Neta,...
          nOfPressureNodesPerElement,NP,material,temperature_el,LS_el,tolerance);

   elseif ismember(ielem,listInt)
       temperature_el = [T_Omega1(T(ielem,:)) T_Omega2(T(ielem,:))];
       LS_el = LS_mesh(T(ielem,:));
       tolerance = InfoProblem.tolerance;
       ngp11 = InfoMesh.ngp_Nits;
       [chigp1,wgp1] = quadrature(InfoMesh.elemType,ngp11);
       [N1,Nxi1,Neta1] = shapeFunctions(InfoMesh.elemType,nOfVelocityNodesPerElement,chigp1);
       [NP1,~,~] = shapeFunctions(elemType,nOfPressureNodesPerElement,chigp1);       
          % element matrix
       [Ke,Ge,fe,muee,tempee,rhoee,Xee] = makeElementMatrix(Xe, ...
          nOfVelocityDofPerElement,chigp1,wgp1,N1,Nxi1,Neta1,...
          nOfPressureNodesPerElement,NP1,material,temperature_el,LS_el,tolerance);
   else
       error('should not go here')
   end

   if plot_up.parameters == 1
       mu = [mu; muee];
       rho = [rho; rhoee];
       temp = [temp; tempee];
       xgp22 = [xgp22; Xee];
   else
       mu = [];
       rho = [];
       temp = [];
       xgp22 = [];
   end 
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

