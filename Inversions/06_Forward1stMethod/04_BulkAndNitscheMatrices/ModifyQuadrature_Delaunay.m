function [chigp,wgp,Area_proportion] = ModifyQuadrature_Delaunay(P1,P2,problem_int,LS_elem,Xe,tolerance)

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
