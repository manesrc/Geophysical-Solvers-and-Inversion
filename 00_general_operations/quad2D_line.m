function [chi_gp,wgp] = quad2D_line(N,A,B)
% receive two points in a 2D-domain and calculate the gauss points and weights 
% using a 1D associated set
% INPUTS
% N : Number of Gauss points as the square of the total amount in each direction
% A : vector of boundaries in isoparametric coordinates (size = [2x1])
% first component initial x value, second component initial y value
% B : vector of boundaries in isoparametric coordinates (size = [2x1])
% first component final x value, second component final y value
% Xe: Nodal coordinates of the element of interest
% OUTPUTS
% chi_gp : matrix with the position of the Gauss points in isoparametric coordinates
% wgp : vector with the weights of the Gauss points

%% code
% isoparametric coordinates
chi_coord = [A(1) B(1)];
eta_coord = [A(2) B(2)];

[chigp,wgp]=quad1D(N,-1,1);
[N,~,~] = shapeFunctions(0,2,chigp); % 0: represent a line element, 2: number of nodes in the elem., chigp position of gauss points in iso
coord1 = N*chi_coord';
coord2 = N*eta_coord';
chi_gp = [coord1 coord2];

return 