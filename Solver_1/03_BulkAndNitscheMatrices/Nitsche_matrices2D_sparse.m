function [G1n,M2s,b1n,m2s,M2s_p] = Nitsche_matrices2D_sparse(InfoMesh, InfoProblem, problem_int,Ar_phys)
% Nitsche_matrices2D_sparse: Assemble Nitsche matrices in 2D using sparse matrices
% Inputs:
%   InfoMesh: structure containing mesh information
%   InfoProblem: structure containing problem information
%   problem_int: 1 - Lithosphere, 2 - Asthenosphere
%   Ar_phys: area fractions of physical domain for each element
% Outputs:
% G1n: Gradient matrix G1n = [ ∫ (∂N/∂x * n) * N d Gamma]
% M2s: Mass matrix M2s =  [ß ∫ N * N d Gamma]
% b1n: Gradient vector b1n = [ ∫ (∂N/∂x * n) * T_LAB d Gamma ]
% m2s: Mass vector m2s = [ ∫ N * T_LAB d Gamma ]
% M2s_p: Mass matrix without considering ß, M2s_p = [ ∫ N * N d Gamma]

%% initialize
% unfold mesh
T = InfoMesh.T;
X = InfoMesh.X;

if problem_int==1
    k = InfoProblem.k1; % diffusivity parameter
    elemNitsche = [InfoMesh.list_cut; InfoMesh.list_edge1]; % elements to calculate Nitsche
else
    k = InfoProblem.k2; % diffusivity parameter
    elemNitsche = [InfoMesh.list_cut; InfoMesh.list_edge2]; % elements to calculate Nitsche
end 

% define num of elements
nelem = size(elemNitsche,1);

% SPARSE
nne = InfoMesh.nne;
% G1
row_G1I = zeros(nne*nne,nelem);
col_G1J = zeros(nne*nne,nelem);
val_G1X = zeros(nne*nne,nelem);
% M2
row_M2I = zeros(nne*nne,nelem);
col_M2J = zeros(nne*nne,nelem);
val_M2X = zeros(nne*nne,nelem);
% M2_unscaled
row_M2uI = zeros(nne*nne,nelem);
col_M2uJ = zeros(nne*nne,nelem);
val_M2uX = zeros(nne*nne,nelem);
% b1
row_b1I = zeros(nne,nelem);
col_b1J = ones(nne,nelem);
val_b1X = zeros(nne,nelem);
% m2
row_m2I = zeros(nne,nelem);
col_m2J = ones(nne,nelem);
val_m2X = zeros(nne,nelem);

%% solve  

% Obtain matrices
for cElem = 1:nelem
  ielem = elemNitsche(cElem,1); % element of interest
  Te = T(ielem,:);  % connect
  Xe = X(Te,:); % coord
  Ar_phys1 = Ar_phys(cElem,2); % find the phys part from list
  assert(ielem == Ar_phys(cElem,1))
  intersection = elemNitsche(cElem,2:5);    % intersction points
  normal2elem = elemNitsche(cElem,end-1:end);   % normal to element
  plot_points_interface = 0;
  [Ge1,Me2,be1,me2,Mee2] = EleMatNitsche(Xe,InfoProblem,InfoMesh,k,problem_int,Ar_phys1,intersection,normal2elem,plot_points_interface);
   % indices
   [mj,mi] = meshgrid(Te,Te);
   % matrix fill 
   % K1
   row_G1I(:,cElem) = mi(:);
   col_G1J(:,cElem) = mj(:);
   val_G1X(:,cElem) = Ge1(:);
   % K2
   row_M2I(:,cElem) = mi(:);
   col_M2J(:,cElem) = mj(:);
   val_M2X(:,cElem) = Me2(:);
   % K2_unscaled
   row_M2uI(:,cElem) = mi(:);
   col_M2uJ(:,cElem) = mj(:);
   val_M2uX(:,cElem) = Mee2(:);
   % vector fill 
   % f1
   row_b1I(:,cElem) = Te(:);
   val_b1X(:,cElem) = be1(:);
  % f2
   row_m2I(:,cElem) = Te(:);
   val_m2X(:,cElem) = me2(:);      
end

ndof_t =size(X,1);

G1n = sparse(row_G1I,col_G1J,val_G1X,ndof_t,ndof_t);
M2s = sparse(row_M2I,col_M2J,val_M2X,ndof_t,ndof_t);
b1n = sparse(row_b1I,col_b1J,val_b1X,ndof_t,1);
m2s = sparse(row_m2I,col_m2J,val_m2X,ndof_t,1);

M2s_p = sparse(row_M2uI,col_M2uJ,val_M2uX,ndof_t,ndof_t);


end