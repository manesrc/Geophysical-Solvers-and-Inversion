function [Mdir,Vdir] = bcFreeSlip(X,nVelocityDof,XP)
% [Mdir,Vdir] = bcFreeSlip(X,nVelocityDof)
% builds the matrix Mdir, and vector Vdir to impose the boundary conditions
% of the Free Slip problem using lagrange multipliers.
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
% put the corners twice
lefNodes = reshape(find(X(:,1) == xmin),[],1);
rigNodes = reshape(find(X(:,1) == xmax),[],1);


% build matrix: [node  bc_value_on_that_node]
% top and bottom: V_y == 0
nodesY = [topNodes; botNodes];
nny = length(nodesY);
bc_Y = [2*nodesY zeros(nny,1)];

% right and left: V_x = 0
nodesX = [lefNodes; rigNodes];
nnx = length(nodesX);
bc_X= [2*nodesX-1 zeros(nnx,1)];


% all boundary conditions together
bc = [bc_X; bc_Y];

% matrix to impose BC using lagrange multipliers
nDirichletBC = size(bc,1);
MdirVel = sparse(nDirichletBC,nVelocityDof);
MdirVel(:,bc(:,1)) = eye(nDirichletBC);
% vector to impose BC using lagrange multipliers
VdirVel = bc(:,2);

%%%%%% pressure boundary condition
npf = find( XP(:,2) == max(XP(:,2)) ) ;      % all the nodes in the Y_max boundary sum zero pressure

% construct a rectangular matrix with ones in npf and condition == 0
nPressDOF = size(XP,1);
PressDir = zeros(1,nPressDOF);
PressDir(npf) = 1;
Press_DirCond = 0;

Mdir = [MdirVel zeros(nDirichletBC,nPressDOF); ...
            zeros(1,nVelocityDof) PressDir];
Mdir = sparse(Mdir);
Vdir = [VdirVel; Press_DirCond];
Vdir = sparse(Vdir);






