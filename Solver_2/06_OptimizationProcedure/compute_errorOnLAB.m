function errorGamma = compute_errorOnLAB(Temp,TLAB,intersected_elements,figurePointsInterface,InfoMesh)
% this code was originally designed for a domain split in two parts Omega1
% and Omega2 above and below the interface respectively.
%
% INPUTS:
% nOfPoints_Elem: the elements of interest are those inters. by the interface and evaluates a 
%                           number of points inside this element (nOfPoints_Elem=2, adds two points 
%                           in each elem inters.)
% DOF_prob: those DOF in the elements inters. by the interface
% problem_int: computes the problem from Omega1 or Omega2 point of view 
%               (normal points out also consider the normal previously computed uses Omega1 point of view)
% figurePointsInterface: to plot the points
% InfoProblem: Material definitions
% InfoMesh: Mesh, connectivities, etc
% PhysicalArea: considers the area of the element inside the physical domain the value is 
%                       normalized to the element area therefore is a ratio (if problem_int == 1, 
%                       is the area of the element in Omega1)
% min_ratio: is a tolerance to consider or not an element to avoid ill-crossing 
%
% OUTPUS: 
% Mat_LABapp:      calculates ∑N_i at points belonging to the interface to project a solution

% CODE: 
% Mesh related information
T = InfoMesh.T;
X = InfoMesh.X;
elemType = InfoMesh.elemType;
nne = InfoMesh.nne;

% points to add between intersection points of interface-element, to evaluate the flux 
nelem2consider = size(intersected_elements,1);
%nOfPoints = nelem2consider; % number of points per element

% Obtain S matrix, the one with the value of the shape functions evaluated on points of the 
% the interface where the jump of flux is measured
%tot_dof = size(X,1);
%matrix_Position = zeros(nOfPoints,tot_dof);

% contador 
errorGamma_elements = 0;
for jj = 1:nelem2consider           % element by element
    % intersection points in isop. coord.
    isop_coord_int1 = intersected_elements(jj,2:3);
    isop_coord_int2 = intersected_elements(jj,4:5);
    iso_coord = [isop_coord_int1; isop_coord_int2];
    % obtain shape functions values
    [N_LABpoints,~,~] = shapeFunctions(elemType,nne,iso_coord);
    % retrieve element connectivities and nodes
    elem_of_interest_jj = intersected_elements(jj,1);
    Te = T(elem_of_interest_jj,:);
    Temp_element = Temp(Te);
    % Temp Gamma:
    TempGamma = N_LABpoints * Temp_element;
    TGamma_mean = sum(TempGamma)/2;
    l_gamma_element = intersected_elements(jj,end); 
    errorGamma_elements = errorGamma_elements + (TGamma_mean - TLAB)^2 * l_gamma_element;
    if figurePointsInterface == 1
        Xe = X(Te,:);
        XPoints = N_LABpoints*Xe;
        figure(120); hold on; 
        scatter(XPoints(:,1),max(X(:,2))-XPoints(:,2),20,'kx')
    end 

end

l_gamma_tot =  sum(intersected_elements(:,end));

admissible_temp_diff = 1*TLAB; 

errorGamma = (1 / (l_gamma_tot * admissible_temp_diff^2)) * errorGamma_elements; 

end 