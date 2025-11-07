function [Mat_LABapp, length_of_influence] = compute_shapeFunctions_atLAB(new_num_points, intersected_elements,figurePointsInterface,InfoMesh)
% computes shape functions at the LAB to measure temperatures. These points are located 
% between the intersection points of the interface and elements and are equally spaced.
% OUTPUT:
% Mat_LABapp: matrix with the shape functions evaluated at the LAB points
% length_of_influence: length of influence associated to each point
% INPUT:
% new_num_points: number of points to add between the intersection points
% intersected_elements: list of intersected elements with the isoparametric coordinates of the intersection points
% figurePointsInterface: flag to plot the points at the interface
% InfoMesh: structure with mesh information


% CODE: 
% Mesh related information
T = InfoMesh.T;
X = InfoMesh.X;
elemType = InfoMesh.elemType;
nne = InfoMesh.nne;

% points to add between intersection points of interface-element, to evaluate the flux 
nelem2consider = size(intersected_elements,1);
nOfPoints = new_num_points*nelem2consider; % number of points per element

% Obtain S matrix, the one with the value of the shape functions evaluated on points of the 
% the interface where the jump of flux is measured
tot_dof = size(X,1);

matrix_Position = zeros(nOfPoints,tot_dof); % one per element and one at each model boundary
length_of_influence = zeros(size(matrix_Position,1),1);

% contador 
matrixRow = 1;
for jj = 1:nelem2consider           % element by element
    % retrieve element connectivities and nodes
    elem_of_interest_jj = intersected_elements(jj,1);
    Te = T(elem_of_interest_jj,:);
    % add points at theboundaries of the model 
    Xe = X(Te,:);    
    % intersection points in isop. coord.
    isop_coord_int1 = intersected_elements(jj,2:3);
    isop_coord_int2 = intersected_elements(jj,4:5);
    % length of influence:
    [N_LABinf,~,~] = shapeFunctions(elemType,nne,[isop_coord_int1; isop_coord_int2]);
    boundaries_XY = N_LABinf*Xe;
    length_interface = sqrt( (boundaries_XY(1,1) - boundaries_XY(2,1))^2 + (boundaries_XY(1,2) - boundaries_XY(2,2))^2 );
    ind_ini = new_num_points*(jj-1)+1;
    ind_fin = new_num_points*jj; 
    length_of_influence(ind_ini:ind_fin) = length_interface/new_num_points; 
    % Calculate the isoparametric position of the midpoint
    xi_new = linspace(isop_coord_int1(1),isop_coord_int2(1),new_num_points+2);
    eta_new = linspace(isop_coord_int1(2),isop_coord_int2(2),new_num_points+2);
    %iso_coord = (isop_coord_int1+isop_coord_int2)/2;
    iso_coord = [xi_new(2:end-1)' eta_new(2:end-1)'];
    % obtain shape functions values
    [N_LABpoints,~,~] = shapeFunctions(elemType,nne,iso_coord);
    % save to the matrix    
    for ii = 1:new_num_points
        matrix_Position(matrixRow,Te) = N_LABpoints(ii,:);
        matrixRow = matrixRow + 1;
    end 
        
    % plot points 
    if figurePointsInterface == 1
        Xe = X(Te,:);
        XPoints = N_LABpoints*Xe;
        figure(120); hold on; 
        scatter(XPoints(:,1),XPoints(:,2),20,'kx')
    end 
end

Mat_LABapp = matrix_Position;   

end 