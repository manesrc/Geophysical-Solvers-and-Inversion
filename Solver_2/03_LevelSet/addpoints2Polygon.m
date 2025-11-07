function densePolygon = addpoints2Polygon(Poly,numPoints)
    %ADDPOINTS2POLYGON adds points in between the vertices of a polygon
    % output:
    %   densePolygon - polygon with added points (M x 2 matrix, M = N*numPoints)
    % input:
    %   Poly - original polygon (Nx2 matrix)
    %   numPoints - number of points to add between each pair of vertices
NumSides = size(Poly,1)-1; % Poly is a closed polygon 
densePolygon = zeros(NumSides*numPoints,2); 

for ii = 1:NumSides
    vertex1 = Poly(ii,:); 
    vertex2 = Poly(ii+1,:);
    points2addX = linspace(vertex1(1),vertex2(1),numPoints+1);
    points2addY = linspace(vertex1(2),vertex2(2),numPoints+1);
    densePolygon((ii-1)*numPoints+1:ii*numPoints,1) = points2addX(1:numPoints);
    densePolygon((ii-1)*numPoints+1:ii*numPoints,2) = points2addY(1:numPoints);
end 

end