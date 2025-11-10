function [chigp,wgp,Area_proportion] = ModifyQuadrature_Delaunay(P1,P2,problem_int,LS_elem,Xe,tolerance)
% This function modifies the quadrature points and weights for an
% element that is cut by the interface, using Delaunay triangulation.
% Inputs:
%   - P1, P2: endpoints of the segment representing the interface within
%             the element (in isoparametric coordinates).
%   - problem_int: Lithosphere is 1 and Asthenosphere is 2.
%   - LS_elem: level set values at the element nodes.
%   - Xe: coordinates of the element nodes.
%   - tolerance: small value to determine the interface location.
%   Outputs:    
%   - chigp: modified quadrature points (in isoparametric coordinates).
%   - wgp: modified quadrature weights.
%   - Area_proportion: proportion of the area of the sub-element.


LSe = LS_elem;

% find X_in and X_out
X_bound = [min(Xe(:,1)) max(Xe(:,1))];
Y_bound = [min(Xe(:,2)) max(Xe(:,2))];
isop_bound = [-1 1];

if problem_int == 1
    logic1 =  LSe > -tolerance;
else
    logic1 = LSe < +tolerance;
end
nodos_relevantes = Xe(logic1,:);

chi_nodes = interp1(X_bound,isop_bound,nodos_relevantes(:,1));
eta_nodes = interp1(Y_bound,isop_bound,nodos_relevantes(:,2));
relevant_nodes_iso = [chi_nodes eta_nodes];

coordenadas = [P1; P2; relevant_nodes_iso];
warning('off')
delauni = delaunayTriangulation(coordenadas);
warning('on')

n_triang = size(delauni.ConnectivityList,1);

elemType_del = 2;       % triangles
num_points = 16;         % in each triangle

[chigp_ref,wgp_ref] = quadrature(elemType_del,num_points);

chigp = zeros(n_triang*num_points,2);
wgp = zeros(1,n_triang*num_points);

contador = 0;
for ii = 1:n_triang
    connectivity_triangle = delauni.ConnectivityList(ii,:);
    coord_triangle = delauni.Points(connectivity_triangle,:);
    [points_real,wgp_real] = mappingTriRefToTri(coord_triangle,chigp_ref,wgp_ref);
    cont_ini = contador*num_points + 1;
    cont_fin = cont_ini + num_points - 1;
    chigp(cont_ini:cont_fin,:) = points_real;
    wgp(cont_ini:cont_fin) = wgp_real';
    contador = contador + 1;
end

Area_proportion = sum(wgp)/4;

end
