    function [Mdir,Vdir] = bcCavity(X,nVelocityDof)
% [Mdir,Vdir] = bcCavity(X,nVelocityDof)
%
% builds the matrix Mdir, and vector Vdir to impose the boundary conditions
% of the cavity problem using lagrange multipliers.
%
% Assumes a square domain.
%

% find domain size
xmin = min(X(:,1)); 
xmax = max(X(:,1));
ymin = min(X(:,2));
ymax = max(X(:,2));

% find nodes on each side
topNodes = reshape(find(X(:,2) == ymax),[],1);
botNodes = reshape(find(X(:,2) == ymin),[],1);
% don't put the corners twice
lefNodes = reshape(find(X(:,1) == xmin & ymin < X(:,2) & X(:,2) < ymax),[],1);
rigNodes = reshape(find(X(:,1) == xmax & ymin < X(:,2) & X(:,2) < ymax),[],1);

% build matrix: [node  bc_value_on_that_node]
nn = length(topNodes);
bctop = [[2*topNodes-1; 2*topNodes] [ones(nn,1); zeros(nn,1)]];

nn = length(botNodes);
bcbot = [[2*botNodes-1; 2*botNodes] zeros(2*nn,1)];

nn = length(lefNodes);
bclef = [[2*lefNodes-1; 2*lefNodes] zeros(2*nn,1)];

nn = length(rigNodes);
bcrig = [[2*rigNodes-1; 2*rigNodes] zeros(2*nn,1)];

% all boundary conditions together
bc = [bctop; bcbot; bclef; bcrig];

% matrix to impose BC using lagrange multipliers
nDirichletBC = size(bc,1);
Mdir = zeros(nDirichletBC,nVelocityDof);
Mdir(:,bc(:,1)) = eye(nDirichletBC);
% vector to impose BC using lagrange multipliers
Vdir = bc(:,2);




