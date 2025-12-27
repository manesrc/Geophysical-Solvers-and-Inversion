function Ge = mgElementMatrix(temp_elem, pres_elem, LS_elem , Xe,velo_element, interface, nnode_el, pespg, N, Nxi, Neta, InfoProblem,InfoMesh, InfoMaterial)

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

LS_gp = N*LS_elem;
temp_gp = N*LS_elem; 
pres_gp = N*pres_elem; 

Ge = ComputeGradMatrix_elemental(temp_gp, pres_gp, LS_gp,Xe,velo_element,N,Nxi,Neta,pespg,nnode_el,cal,InfoMaterial);

end


function Ge = ComputeGradMatrix_elemental(tempgp, presgp, LSgp,Xe,velocity,N,Nxi,Neta,pespg,numberOfNodes,cal,material)

numberOfGaussPoints = length(pespg); 
Ge = zeros(numberOfNodes, numberOfNodes); 

rho = DensityAtGaussPoints(LSgp,tempgp,presgp,293,material);
%rho = 3300 * ones(size(rho));

    for igaus = 1:numberOfGaussPoints 
        if pespg(igaus) ~= 0
            N_igaus = N(igaus,:);
            %X_igaus = N_igaus*Xe; 
            velo_gp = N_igaus * velocity;
            jacob = [Nxi(igaus,:)*Xe(:,1)  Nxi(igaus,:)*Xe(:,2)
                Neta(igaus,:)*Xe(:,1) Neta(igaus,:)*Xe(:,2)]; 
           dvolu = pespg(igaus) * det(jacob); 
           res = jacob\[Nxi(igaus,:); Neta(igaus,:)]; 
           u_gradT = res' * velo_gp';
           %
           Ge = Ge + rho(igaus) * cal * N_igaus' * u_gradT' *dvolu; 
        end 
    end

end