function [pts,w] = mappingTriRefToTri(verticesXv,gaussPoints,gaussWeights)
%MAPPINGTRIREFTOTRI Maps Gauss points and weights from reference triangle to an
%   arbitrary triangle defined by its vertices (for Delaunay triangulation).
%   Inputs:
%       verticesXv - 3x2 matrix where each row represents the (x,y
%                     coordinates of a vertex of the triangle.
%       gaussPoints - Nx2 matrix of Gauss points in the reference triangle.
%       gaussWeights - Nx1 vector of Gauss weights in the reference triangle.
%   Outputs:
%       pts - Nx2 matrix of mapped Gauss points in the arbitrary triangle.
%       w - Nx1 vector of mapped Gauss weights in the arbitrary triangle.

% From (-1,-1) - (1,-1) - (-1,1) 

a = verticesXv(1,:);
b = verticesXv(2,:);
c = verticesXv(3,:);

%pts = (gaussPoints + 1 )*0.5 ; 
pts = gaussPoints;
%w = gaussWeights/4.0 ;
A = [  a(1) - c(1) , b(1) - c(1) ; a(2) - c(2) , b(2) - c(2) ] ; 

nOfPts = size(pts,1);
for i = 1 : nOfPts
    pts(i,:) = (A* (pts(i,:)'))' + c;   
end


TriArea = polyarea(verticesXv(:,1),verticesXv(:,2));
w = gaussWeights/sum(gaussWeights)*TriArea ;


end