function [Ke,fe,Ar_phys_el] = mkElementMatrix(Xe, numberOfNodes,pespg, N, Nxi, Neta, ...
    InfoProblem,problem_int,interface,InfoMesh)
%% elemental matrix calculation
%Xe, nnode, pospg, pespg, N, Nxi, Neta, InfoProblem,problem_int,interface,InfoLAB


if problem_int == 1
   conductivity = InfoProblem.k1;
   source = InfoProblem.s1;
else 
   conductivity = InfoProblem.k2;
   source = InfoProblem.s2;
end 

% measure the area of the physical part of the elemeTnt crossed by LAB
if interface.elem ==  1
    % Use Level set (coord. [X,Y]) and compute gauss points, weight and
    % parts of the element in the domain [fact_phys]
    LS_elem = interface.LS_elem; 
    if interface.tryDelaunay == 1
        P1 = interface.P_in_xieta; 
        P2 = interface.P_out_xieta;
        [chigp1,wgp1,Ar_phys_el] = ModifyQuadrature_Delaunay(P1,P2,problem_int,LS_elem,Xe,InfoProblem.tolerance);
    else
        [chigp1,wgp1,Ar_phys_el] = ModifyQuadrature_LevelSet(InfoMesh,problem_int,LS_elem,InfoProblem.tolerance);
    end 

    % calculate shape functions and derivatives 
    [N_mod,Nxi_mod,Neta_mod] = shapeFunctions(InfoMesh.elemType,InfoMesh.nne,chigp1);
    
    if interface.plot == 1
        % plot points to illustrate the integration
        points_considered = N_mod*Xe; 
        if problem_int == 1
            figure(120); hold on; scatter(points_considered(wgp1~=0,1),points_considered(wgp1~=0,2),20,'black','^','filled')
        else
            figure(120); hold on; scatter(points_considered(wgp1~=0,1),points_considered(wgp1~=0,2),20,'red','o','filled')
        end 
    end
    % calculate elemental matrix
    [Ke,fe] = Compute_element_matrix(Xe,N_mod,Nxi_mod,Neta_mod,wgp1,numberOfNodes,conductivity,source);
else 
    [Ke,fe] = Compute_element_matrix(Xe,N,Nxi,Neta,pespg,numberOfNodes,conductivity,source);
    Ar_phys_el = []; % in this element fact_phys is irrelevant
end 
end 

function [Ke,fe] = Compute_element_matrix(Xe,N,Nxi,Neta,pespg,numberOfNodes,conductivity,source)
numberOfGaussPoints = length(pespg); 
Ke = zeros(numberOfNodes, numberOfNodes); 
fe = zeros(numberOfNodes, 1);
    for igaus = 1:numberOfGaussPoints 
        if pespg(igaus) ~= 0
            jacob = [Nxi(igaus,:)*Xe(:,1)  Nxi(igaus,:)*Xe(:,2)
                Neta(igaus,:)*Xe(:,1) Neta(igaus,:)*Xe(:,2)]; 
           dvolu = pespg(igaus) * det(jacob); 
           res = jacob\[Nxi(igaus,:); Neta(igaus,:)]; 
           Nx = res(1,:); 
           Ny = res(2,:);
           %
           Ke = Ke + conductivity*(Nx'*Nx + Ny'*Ny)*dvolu; 
           fe = fe + source * N(igaus,:)' *dvolu; 
        end 
    end
end
