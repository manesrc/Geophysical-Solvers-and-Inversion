function [G_matrix,matrixS] = compute_Grad_interphase(nOfPoints_Elem,DOF_prob,problem_int,figurePointsInterface,InfoProblem,InfoMesh,Physical_area,min_ratio)
% compute_Grad_interphase: computes the gradient matrix at the interphase
% between two subdomains Lithosphere or Asthenosphere (1 or 2)
%
% INPUT:   
% nOfPoints_Elem        : number of points to evaluate the flux per element
% DOF_prob              : degrees of freedom of the problem to consider
% problem_int          : 1 or 2, depending on the subdomain to consider
% figurePointsInterface : 1 to plot the points where the flux is evaluated
% InfoProblem          : structure with information of the problem
% InfoMesh             : structure with information of the mesh
% Physical_area        : matrix with the area of each element in each subdomain
% min_ratio            : minimum ratio of area to consider an element well supported
%
% OUTPUT:
% G_matrix             : gradient matrix evaluated at points at the interphase
% matrixS              : shape functions matrix evaluated at points at the interphase

% Mesh related information
T = InfoMesh.T;
X = InfoMesh.X;
elemType = InfoMesh.elemType;
nne = InfoMesh.nne;


if isempty(InfoMesh.list_cut) ~= 1
    % this list contains the elements that have small areas in Omega_1 and/or
    % Omega_2
    elements_illSupport = [];
    elements_in_list = 1:1:length(InfoMesh.list_cut);
    for ii = 1:size(InfoMesh.list_cut,1)
        elem_interest = InfoMesh.list_cut(ii,1);
        logical1 = Physical_area(:,1)==elem_interest;
        Elem_phys_area = Physical_area(logical1,2);
        if Elem_phys_area < min_ratio || Elem_phys_area > (1-min_ratio)
            elements_illSupport = [elements_illSupport; elements_in_list(logical1)];
        end 
    end 
    list_cut_trim = InfoMesh.list_cut;
    list_cut_trim(elements_illSupport,:) = [];
else
    list_cut_trim = [];
end 

if problem_int == 1
    % define the points to evaluate the flux
    %intersected_elements = [InfoMesh.list_cut; InfoMesh.list_edge1]; % elements and intersections with interface
    intersected_elements = [list_cut_trim; InfoMesh.list_edge1]; % elements and intersections with interface
    k_diffussivity = InfoProblem.k1;
    normal_sign = 1;
elseif problem_int == 2
    %intersected_elements = [InfoMesh.list_cut; InfoMesh.list_edge2]; % elements and intersections with interface
    intersected_elements = [list_cut_trim; InfoMesh.list_edge2]; % elements and intersections with interface
    k_diffussivity = InfoProblem.k2;
    normal_sign = -1;
end 


% points to add between intersection points of interface-element, to evaluate the flux 
nelem2consider = size(intersected_elements,1);
nOfPointsFlux = nOfPoints_Elem*nelem2consider; % number of points per element

% Obtain S matrix, the one with the value of the shape functions evaluated on points of the 
% the interface where the jump of flux is measured
tot_dof = size(InfoMesh.X,1);
matrix_Gradiente = zeros(nOfPointsFlux,tot_dof);
matrix_Position = zeros(nOfPointsFlux,tot_dof);

% contador 
matrixRow = 0;
for jj = 1:nelem2consider           % element by element
    % intersection points in isop. coord.
    isop_coord_int1 = intersected_elements(jj,2:3);
    isop_coord_int2 = intersected_elements(jj,4:5);
    % normal in the element of interest
    normal_interest = normal_sign*intersected_elements(jj,end-1:end);
    % add points between intersection points
    mid_points_chi = linspace(isop_coord_int1(1),isop_coord_int2(1),nOfPoints_Elem+2);
    mid_points_eta = linspace(isop_coord_int1(2),isop_coord_int2(2),nOfPoints_Elem+2);
    iso_coord = [mid_points_chi(2:end-1)' mid_points_eta(2:end-1)'];
    % obtain shape functions values
    [N_LABpoints,Nxi_LABp,Neta_LABp] = shapeFunctions(elemType,nne,iso_coord);

    % retrieve element connectivities and nodes
    elem_of_interest_jj = intersected_elements(jj,1);
    Te = T(elem_of_interest_jj,:);
    Xe = X(Te,:); 
    % save to the matrix
    for qq=1:nOfPoints_Elem
        mat1 = [Nxi_LABp(qq,:); Neta_LABp(qq,:)];
        Jacob = mat1*Xe; 
        res = Jacob\mat1;
        Nx_qq = res(1,:);
        Ny_qq = res(2,:);
        Contracted_Gradient = [Nx_qq' Ny_qq']*normal_interest'; 
        matrix_Gradiente(matrixRow+qq,Te) = Contracted_Gradient;  % put the values of the SF in each row
        matrix_Position(matrixRow+qq,Te) = N_LABpoints(qq,:);
    end 
    matrixRow = matrixRow + nOfPoints_Elem;
    if figurePointsInterface == 1
        Xe = InfoMesh.X(Te,:);
        XPoints = N_LABpoints*Xe;
        figure(120); hold on; 
        scatter(XPoints(:,1),XPoints(:,2),20,'kx')
    end 

end

G_matrix = k_diffussivity*matrix_Gradiente(:,DOF_prob);
%G_matrix = matrix_Gradiente;
matrixS = matrix_Position; %(:,DOF_prob);
%matrixS1 =matrix_Position(:,DOF_prob); 
assert(matrixRow == size(G_matrix,1))
assert(nelem2consider * nOfPoints_Elem == matrixRow)

end 