function [Mdir,Vdir] = bcFreeSlip(X,nVelocityDof,XP)
% Applies free slip boundary conditions using Lagrange multipliers. Assumes a square domain.
% outputs:
% Mdir : matrix to impose the BCs using Lagrange multipliers    
% Vdir : vector to impose the BCs using Lagrange multipliers
% inputs:
% X : coordinates of the velocity nodes
% nVelocityDof : number of velocity degrees of freedom
% XP : coordinates of the pressure nodes


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
bc_Y = [2*nodesY zeros(nny,1)];     % 2*nodesY (means the Y component of the velocity of that node)

% right and left: V_x = 0
nodesX = [lefNodes; rigNodes];
nnx = length(nodesX);
bc_X= [2*nodesX-1 zeros(nnx,1)]; % 2*nodesX-1 (means the X component of the velocity of that node)


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

end






