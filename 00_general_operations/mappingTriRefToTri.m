function [pts,w] = mappingTriRefToTri(verticesXv,gaussPoints,gaussWeights)
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