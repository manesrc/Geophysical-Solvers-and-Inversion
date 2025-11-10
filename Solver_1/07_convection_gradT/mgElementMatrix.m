function Ge = mgElementMatrix(Xe,velo_element, interface, nnode_el, pespg, N, Nxi, Neta, InfoProblem,InfoMesh, InfoMaterial)
% mgElementMatrix: computes the elemental matrix for the convection term with gradT
% INPUT:
%   Xe: coordinates of the element nodes
%   velo_element: velocity at the element nodes
%   interface: structure with information about the interface
%   nnode_el: number of nodes in the element
%   pespg: weights of the Gauss points
%   N, Nxi, Neta: shape functions and their derivatives at Gauss points
%   InfoProblem: structure with problem information
%   InfoMesh: structure with mesh information
%   InfoMaterial: structure with material properties
% OUTPUT:
%   Ge: elemental matrix for the convection term with gradT 

rho = 1;        % adimensional calculation
cal = InfoMaterial.calorific; % adimensional calculation

if interface.elem == 1
    % compute the integration in an element cut by the interface
    % Use Level set (coord. [X,Y]) and compute gauss points, weight and parts of 
    % the element in the domain [fact_phys]
    LS_elem = interface.LS_elem; 
    problem_int = 2;
    if interface.tryDelaunay == 1
        P1 = interface.P_in_xieta; 
        P2 = interface.P_out_xieta;
        [chigp1,pespg,~] = ModifyQuadrature_Delaunay(P1,P2,problem_int,LS_elem,Xe,InfoProblem.tolerance);
    else
        [chigp1,pespg,~] = ModifyQuadrature_LevelSet(InfoMesh,problem_int,LS_elem,InfoProblem.tolerance);
    end 
    % calculate shape functions and derivatives 
    [N,Nxi,Neta] = shapeFunctions(InfoMesh.elemType,InfoMesh.nne,chigp1);
end 

Ge = ComputeGradMatrix_elemental(Xe,velo_element,N,Nxi,Neta,pespg,nnode_el,rho,cal);

end


function Ge = ComputeGradMatrix_elemental(Xe,velocity,N,Nxi,Neta,pespg,numberOfNodes,rho,cal)

numberOfGaussPoints = length(pespg); 
Ge = zeros(numberOfNodes, numberOfNodes); 

    for igaus = 1:numberOfGaussPoints 
        if pespg(igaus) ~= 0
            N_igaus = N(igaus,:);
            velo_gp = N_igaus * velocity;
            jacob = [Nxi(igaus,:)*Xe(:,1)  Nxi(igaus,:)*Xe(:,2)
                Neta(igaus,:)*Xe(:,1) Neta(igaus,:)*Xe(:,2)]; 
           dvolu = pespg(igaus) * det(jacob); 
           res = jacob\[Nxi(igaus,:); Neta(igaus,:)]; 
           u_gradT = res' * velo_gp';
           %
           Ge = Ge + rho * cal * N_igaus' * u_gradT' *dvolu; 
        end 
    end

end