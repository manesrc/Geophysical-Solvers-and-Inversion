function [list_Omega1,list_Omega2,list_cut_elem,list_edge_elem1,list_edge_elem2, LABpoints] = findElemCrossedLS(tol,plot1,InfoMesh,InfoLAB,InfoProblem)
% This function classifies the elements of the mesh according to the level set function
% Outputs:
% list_Omega1: list of elements in Omega1
% list_Omega2: list of elements in Omega2
% list_cut_elem: list of elements cut by the LAB + intersection points + normal to Omega1
% list_edge_elem1: list of ill-crossed elements in Omega1 + intersection points + normal to Omega1
% list_edge_elem2: list of ill-crossed elements in Omega2 + intersection points + normal to Omega1
% LABpoints: points in the LAB to evaluate magnitudes at the LAB
% Inputs:
% tol: tolerance to consider level set values as zero
% plot1: flag to plot the crossed elements (1: yes, 0: no)
% InfoMesh: structure with mesh information
% InfoLAB: structure with LAB information
% InfoProblem: structure with problem information   

% level set
LS = InfoMesh.LS_new;
% mesh matrices
X = InfoMesh.X;
T = InfoMesh.T;

% data to save
list_Omega1 = [];
list_Omega2 = [];
list_cut_elem = [];
list_edge_elem1 = [];
list_edge_elem2 = [];

processed_elements = false(size(T,1),1);


for jj = 1:size(T,1)
    % Skip if the element has already been processed
    if processed_elements(jj)
        continue;
    end

    Te = T(jj,1:4);     % simplified element [linear]
    LSe = LS(Te,:);
    % apply tolerances
    for ii = 1:size(LSe,1)
        if abs(LSe(ii)) < tol
            LSe(ii) = 0;
        end 
    end 
    % 
    count_pos = sum(LSe > 0);
    count_neg = sum(LSe < 0);
    count_zeros = sum(LSe == 0);

    count_tot = length(LSe);
    
    assert((count_pos+count_neg+count_zeros) == count_tot)
    
    % elements in Omega2 
    if count_neg == count_tot || (count_neg == count_tot-1 && count_zeros == 1)
        list_Omega2 = [list_Omega2; jj];      % element in Omega 2
    
    % elements in Omega1
    elseif count_pos == count_tot || (count_pos == count_tot-1 && count_zeros == 1)
        list_Omega1 = [list_Omega1; jj];      % element in Omega 1
    
    % elements crossed by LAB 
    elseif count_pos > 0 && count_neg > 0
        % intersection points
        Xe_corners = X(Te,:);   % X,Y points of the corners of the element
        [intersec_1,intersec_2] = findElemInter_intersection(Xe_corners,LSe);   % calculates the points (isop. coord.) of the intersection element-interface
        % normal pointing out of Omega1 considering the intersection points (interface linearly approx. in the element)
        normal2Omega1 = compute_normal_Omega1(intersec_1,intersec_2,Xe_corners,LSe);
        % list of element crossed by interface + geometric characteristics of interest
        list_cut_elem = [list_cut_elem; jj intersec_1 intersec_2 normal2Omega1];      
    
    % elements with an arista coinciding with LAB (ill-crossed)
    elseif (count_zeros == 2) && (count_pos == 2 || count_neg == 2)      % case line coinciding with arista intersection points (coincident with zeros)
        % Xe_corners = X(Te,:);   % X,Y points of the corners of the element
        % [intersec_1,intersec_2] = findElemInter_intersection(Xe_corners,LSe);   % calculates the points (isop. coord.) of the intersection element-interface
        % check_intersection = [intersec_1 intersec_2];
        % if any(abs(check_intersection)~=1)
        %     error('the intersection points are not aligned with the element boundaries')
        % end 
        % % normal pointing out of Omega1 considering the intersection points (interface linearly approx. in the element)
        % normal2Omega1 = compute_normal_Omega1(intersec_1,intersec_2,Xe_corners,LSe);
        % 
        % 
        % % if ( (abs(normal2Omega1(1)) == 0 || abs(normal2Omega1(2)) == 0) && (abs(normal2Omega1(1)) == 1 || abs(normal2Omega1(2)) == 1) ) ~= 1
        % %     error('the element normal do not match with the intersection')
        % % end 

        % Case: line coinciding with arista (only valid for [1 2], [2 3], [3 4], [4 1])
        % Find which nodes are zero
        zero_nodes = find(LSe == 0);
        
        % Check if these nodes form a valid arista
        valid_aristas = [1 2; 2 3; 3 4; 4 1];
        if any(ismember(sort(zero_nodes)', valid_aristas, 'rows'))
            % This is a valid arista case
            % Process the element as usual (existing code here)
            
            Xe_corners = X(Te,:);   % X,Y points of the corners of the element
            [intersec_1,intersec_2] = findElemInter_intersection(Xe_corners,LSe); % calculates the points (isop. coord.) of the intersection element-interface
            
            % Normal pointing out of Omega1 considering the intersection points (interface linearly approx. in the element)
            normal2Omega1 = compute_normal_Omega1(intersec_1,intersec_2,Xe_corners,LSe);
            
            % Classify to appropiate list
            if count_neg == 2
                list_edge_elem2 = [list_edge_elem2; jj intersec_1 intersec_2 normal2Omega1];
                % Find the element in contact with this one, belonging to Omega1
                elem_in_domain = 2;
                [other_elem, intersection1, intersection2, normal2Omega1] = find_elem_in_contact(jj, elem_in_domain, normal2Omega1, InfoMesh, list_Omega1, list_Omega2, list_cut_elem);
                list_edge_elem1 = [list_edge_elem1; other_elem intersection1 intersection2 normal2Omega1];
                % Mark the neighboring element as processed
                processed_elements(other_elem) = true;
            else
                list_edge_elem1 = [list_edge_elem1; jj intersec_1 intersec_2 normal2Omega1];
                % Find the element in contact with this one, belonging to Omega2
                elem_in_domain = 1;
                [other_elem, intersection1, intersection2, normal2Omega1] = find_elem_in_contact(jj, elem_in_domain, normal2Omega1, InfoMesh, list_Omega1, list_Omega2, list_cut_elem);
                list_edge_elem2 = [list_edge_elem2; other_elem intersection1 intersection2 normal2Omega1];
                % Mark the neighboring element as processed
                processed_elements(other_elem) = true;
            end
        elseif count_neg == 2
            list_Omega2 = [list_Omega2; jj];      % element in Omega 2
        elseif count_pos == 2
            list_Omega1 = [list_Omega1; jj];      % element in Omega 1
        end
            % if count_neg == 2
            %     list_edge_elem2 = [list_edge_elem2; jj intersec_1 intersec_2 normal2Omega1];
            %     % find element in contact with the previous one, belonging to Omega1
            %     elem_in_domain = 2;
            %     [other_elem, intersection1, intersection2, normal2Omega1] = find_elem_in_contact(jj, elem_in_domain,normal2Omega1,InfoMesh,list_Omega1,list_Omega2,list_cut_elem);
            %     list_edge_elem1 = [list_edge_elem1; other_elem intersection1 intersection2 normal2Omega1];
            %     % Mark the neighboring element as processed
            %     processed_elements(other_elem) = true;
            % else
            %     list_edge_elem1 = [list_edge_elem1; jj intersec_1 intersec_2 normal2Omega1];
            %     % find element in contact with the previous one, belonging to Omega2
            %     elem_in_domain = 1;
            %     [other_elem, intersection1, intersection2, normal2Omega1] = find_elem_in_contact(jj, elem_in_domain,normal2Omega1,InfoMesh,list_Omega1,list_Omega2,list_cut_elem);
            %     list_edge_elem2 = [list_edge_elem2; other_elem intersection1 intersection2 normal2Omega1];
            %     % Mark the neighboring element as processed
            %     processed_elements(other_elem) = true;
            % end 

    else
        error('Si entra acá: FALTAN COSAS')
    end 
end 

if size(list_edge_elem2,1) ~= size(list_edge_elem1,1)
    disp('elements in edge domain 2')
    list_edge_elem2
    disp('elements in edge domain 1')
    list_edge_elem1 
    error('problem with elements')
end 

%% assert the all elements are classified just once
num_of_elements_classified = ( length(list_Omega1) + length(list_Omega2) + size(list_cut_elem,1) + size(list_edge_elem1,1) + size(list_edge_elem2,1) ) ;

if num_of_elements_classified ~= size(T,1)
    keyboard
    % check if an element has been classificated in two lists at the same time:
    % Omega1 and Omega2
    list_Omega1 = check_duplicated_classification(list_Omega1,list_Omega2);
    % well-cut and Omega1
    list_cut_elem = check_duplicated_classification(list_cut_elem,list_Omega1);
    % well-cut and Omega2
    list_cut_elem = check_duplicated_classification(list_cut_elem,list_Omega2);
    % ill-cut1 and well-cut
    list_edge_elem1 = check_duplicated_classification(list_edge_elem1,list_cut_elem);
    % ill-cut2 and well-cut
    list_edge_elem2 = check_duplicated_classification(list_edge_elem2,list_cut_elem);
    % ill-cut1 and ill-cut2
    list_edge_elem1 = check_duplicated_classification(list_edge_elem1,list_edge_elem2);
end 

%% save

list_of_elements = [list_cut_elem; list_edge_elem1];
plot_elem = 0; 
LABpoints = plot_interphase_from_intersection(list_of_elements,InfoMesh.elemType,InfoMesh.nne,X,T,plot_elem);

%% plot crossed elements 
if plot1 == 1
    plot_crossed_elements(list_cut_elem,list_edge_elem1,list_edge_elem2,LABpoints,InfoMesh,InfoProblem)
end 

end

function [chieta1, chieta2] = findElemInter_intersection(Xe, LSe)
    ind = 1:1:length(LSe);
    % determine points of the interface entering and leaving the element
    % pos. and neg. LS nodes 
    LS_pos = LSe(LSe > 0);
    LS_neg = LSe(LSe < 0);
    % if any LS == 0 -> there's an intersection point 
    intersPoint_LS0 = ind(LSe == 0);


    if length(intersPoint_LS0) == 2
        % int. points in X,Y coord.
        x_int1 = Xe(ind(intersPoint_LS0(1)),:);
        x_int2 = Xe(ind(intersPoint_LS0(2)),:);
        % find isop. coord. 
        [chi1,eta1] = isoparametric_coordinates(Xe,x_int1(1),x_int1(2));
        [chi2,eta2] = isoparametric_coordinates(Xe,x_int2(1),x_int2(2));
        % points of interest
        chieta1 = [chi1 eta1];
        chieta2 = [chi2 eta2];
    elseif length(intersPoint_LS0) == 1
        % int. point in X,Y coord.
        x_int1 = Xe(ind(intersPoint_LS0),:);
        % find isop. coord. 
        [chi1,eta1] = isoparametric_coordinates(Xe,x_int1(1),x_int1(2));
        % point of interest 1
        chieta1 = [chi1 eta1];
        % second point: [find the node with different sign LS]
        indie = ind; 
        if length(LS_pos) == 1
            indi_exist = ind(intersPoint_LS0);
            indi = ind(LSe > 0);
            indie(indie == indi_exist) = []; indie(indie==indi) = []; 
            if ismember(indi+1,indie)
                indi1 = indi+1;
                x_int2 = Xe(indi1,:);
            elseif ismember(indi-1,indie)
                indi1 = indi-1;
                x_int2 = Xe(indi1,:);
            else
                bounds1 = [1 4];
                indi1 = bounds1(ismember(bounds1,indie));
                x_int2 = Xe(indi1,:);
            end 
            x_arista1 = [Xe(indi,:); x_int2];
            LS1 = [LSe(indi); LSe(indi1)];
            chieta2 = point_in_arista(x_arista1,LS1,Xe);
        else
            indi_exist = ind(intersPoint_LS0);
            indi = ind(LSe<0);
            indie(indie == indi_exist) = []; indie(indie==indi) = []; 
            if ismember(indi+1,indie)
                indi1 = indi+1;
                x_int2 = Xe(indi1,:);
            elseif ismember(indi-1,indie)
                indi1 = indi-1;
                x_int2 = Xe(indi1,:);
            else
                bounds1 = [1 4];
                indi1 = bounds1(ismember(bounds1,indie));
                x_int2 = Xe(indi1,:);
            end 
            x_arista1 = [Xe(indi,:); x_int2];
            LS1 = [LSe(indi); LSe(indi1)];
            chieta2 = point_in_arista(x_arista1,LS1,Xe);
        end
    else
        if (length(LS_pos) == 1) || (length(LS_neg) == 1)                                                                 
            % which point is with a different sign 
            if length(LS_pos) == 1
                indi = ind(LSe > 0);
            else
                indi = ind(LSe < 0);
            end 
            % other two points sharing an arista
            indi1 = indi-1; 
            indi2 = indi+1;
            if indi1 == 0
                indi1 = 4;
            elseif indi2 == 5
                indi2 = 1;
            end 
            % first intersection point: arista 1
            x_arista1 = [Xe(indi,:); Xe(indi1,:)]; % points in arista 1
            LS1 = [LSe(indi); LSe(indi1)];                     % values of the level set in the nodes of the arista 
            chieta1 = point_in_arista(x_arista1,LS1,Xe);     % intersection interface - arista 1  
            % second intersection point: arista 2
            x_arista2 = [Xe(indi,:); Xe(indi2,:)]; % points in arista 2
            LS2 = [LSe(indi); LSe(indi2)];                     % values of the level set in the nodes of the arista
            chieta2 = point_in_arista(x_arista2,LS2,Xe);     % intersection interface - arista 2
        else
            % element crossed from left-right or up-down
            x11 = Xe(LSe > 0,:); % this matrix has two points on the same side of the LS
            if x11(1,1) == x11(2,1)     % the LS is positive/negative on the left/right sides
                x_arista1 = [Xe(1,:); Xe(2,:)];
                LS1 = [LSe(1) LSe(2)];
                x_arista2 = [Xe(3,:); Xe(4,:)];
                LS2 = [LSe(3) LSe(4)];
            else        % the LS is positive/negative on the top/bottom sides
                x_arista1 = [Xe(1,:); Xe(4,:)];
                LS1 = [LSe(1); LSe(4)];
                x_arista2 = [Xe(2,:); Xe(3,:)];
                LS2 = [LSe(2); LSe(3)];
            end
            chieta1 = point_in_arista(x_arista1,LS1,Xe);     % intersection interface - arista 1  
            chieta2 = point_in_arista(x_arista2,LS2,Xe);     % intersection interface - arista 2
    
            %x12 = Xe_simpl(LSe_simpl < 0,:); % MAL
            %error('FALTA TERMINAR PAPAAA!!!!!!!!!')
            %LS1 = LSe
        end 
    end 
end
function chieta = point_in_arista(x1,LS1,X1)
    % OUTPUT: 
    % chieta = [chi eta] of intersection interface - arista
    % INPUT:
    % x1 = two points defining the arista
    % LS1 = the values of the level set in the arista
    % X1 = four points of the corners of the element
    
    if LS1(1) == 0 || LS1(2) == 0
        if LS1(1) == 0
            % condicion 1 
            [chi1,eta1] = isoparametric_coordinates(X1,x1(1,1),x1(1,2));
            chieta = [chi1 eta1];
        else 
            [chi1,eta1] = isoparametric_coordinates(X1,x1(2,1),x1(2,2));
            chieta = [chi1 eta1];
        end 
    else
        if x1(1,1) == x1(2,1)       % arista is parallel to Y-axis
            vect1 = [x1(1,2) x1(2,2)];
            p = polyfit(vect1,LS1,1);
            Ycoord_LS0 = roots(p);
            % obtain the isoparametric position 
            [chi1,eta1] = isoparametric_coordinates(X1,x1(1,1),Ycoord_LS0);
            chieta = [chi1 eta1];
        else        % arista is parallel to X-axis
            vect1 = [x1(1,1) x1(2,1)];
            p = polyfit(vect1,LS1,1);
            Xcoord_LS0 = roots(p);
            % obtain the isoparametric position 
            [chi1,eta1] = isoparametric_coordinates(X1,Xcoord_LS0,x1(1,2));
            chieta = [chi1 eta1];
        end 
    end 
    
end 
function normal2Omega1_element = compute_normal_Omega1(point1,point2,Xe,LSe)
    % output: normal to domain Omega_1
    % inputs: point1 and point2 interface goes into the elem. point1
    %                                           interface comes out of the elem, point2 
    %             Xe: points of the corners of the element
    %             LSe: values of the level set in the corners of the element
    % from the isoparametric coordinates, compute the X,Y coord.
    X_bound = [min(Xe(:,1)) max(Xe(:,1))];      
    Y_bound = [min(Xe(:,2)) max(Xe(:,2))];
    isop_bounds = [-1 1];
    chi_data = [point1(1) point2(1)];
    x_data = interp1(isop_bounds,X_bound,chi_data);
    eta_data = [point1(2) point2(2)];
    y_data = interp1(isop_bounds,Y_bound,eta_data);
    % intersection points in X,Y coord.
    point1_xy = [x_data(1) y_data(1)];
    point2_xy = [x_data(2) y_data(2)];
    % interface line
    tramo_12 = (point2_xy - point1_xy);         % Point_out - Point_in 
    % normal to the interface
    normal2Omega = [tramo_12(2) -tramo_12(1)]/norm(tramo_12);
    
    %% now check the normal2Omega is contrary to the gradient of the Level Set
    % point in the middle of interface line
    point_cent_interf_chieta = (point1+point2)/2;
    % shape functions to compute the gradient
    elemType = 1; % quadrilaterals
    nElemNodes = size(Xe,1);    % linear quads
    
    [~,Nxi,Neta] = shapeFunctions(elemType,nElemNodes,point_cent_interf_chieta);
    gradN_iso = [Nxi;Neta];
    Jacob = gradN_iso * Xe;
    grad_LS = (Jacob\gradN_iso) * LSe ;
    
    if (normal2Omega * grad_LS) < 0 
        normal2Omega1_element = normal2Omega;
    else
        normal2Omega1_element = -normal2Omega;
    end 
end 
function addElem(X,arg)
   %ix = [1 2 3 4 1 5 6 7 8 5 6 2 3 7 8 4]; % malla lineal
   ix = [1 2 3 4 1];
   plot(X(ix,1),X(ix,2),arg{:})
end
function addElem2(X,arg)
   ix = [1 2 3 4 1]; % malla lineal
%    ix = [1 2 3 4 1 13 14 15 16 13 14 2 3 15 16 4];
   plot(X(ix,1),X(ix,2),arg{:})
end
function [chi,eta] = isoparametric_coordinates(Xe,x,y)
    % calculate the position of the x,y coordinates in isoparametric
    % description inside the element with Xe coordinates
    Xboundaries = [min(Xe(:,1)) max(Xe(:,1))];
    Yboundaries = [min(Xe(:,2)) max(Xe(:,2))];
    % obtain isoparametric position of the intersections 
    iso_coord = [-1 1]; 
    % calculate isoparametric coord. 
    chi = interp1(Xboundaries, iso_coord,x);
    eta = interp1(Yboundaries, iso_coord,y);
end
function LABpoints2 = plot_interphase_from_intersection(list_of_cut_elements,elemType,nne,X,T,plot2)
    num_elements = size(list_of_cut_elements,1);
    ymax = max(X(:,2));
    LABpoints2 = zeros(2,2*num_elements);
    
    for ii = 1:num_elements
        elem_ii = list_of_cut_elements(ii,1);
        Xe = X(T(elem_ii,:),:);
        xigp1 = list_of_cut_elements(ii,2:3);
        xigp2 = list_of_cut_elements(ii,4:5);
        [Nii,~,~] = shapeFunctions(elemType,nne,[xigp1; xigp2]);
        Xgpii = Nii*Xe;
        Xgpii(:,2) = ymax - Xgpii(:,2);
        ini_index = 2*(ii-1)+1;
        fin_index = ini_index+1;
        LABpoints2(:,ini_index:fin_index) = Xgpii'; 
        if plot2 == 1
            hold on
            plot(Xgpii(:,1),ymax-Xgpii(:,2),'Color',[0.8500 0.3250 0.0980],'Marker','s','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize',1);
        end 
    end 
    plotSortedPoints = 0; 
    [LABx,LABy] = sortLABPoints2plot(LABpoints2(1,:),LABpoints2(2,:),101,plotSortedPoints);
    LABpoints2 = [LABx'; LABy']; 
end 
function [other_elem, intersection1, intersection2, normal2Omega1,list1_new,list2_new,list_cut_new] = find_elem_in_contact(elem_num, belong2domain,normal2Omega1_original,InfoMesh,list1,list2,list_cut) %,list_Omega1,list_Omega2,list_cut_elem
    % this function aims at determining the element in contact with the
    % interface to the other side to the one that has been already found
    X = InfoMesh.X;
    T = InfoMesh.T;
    nel_x = InfoMesh.nel_x; 
    if abs(normal2Omega1_original(1)) == 1   
        % arista coincides with lateral element boundary
        % element in Omega2
        if ((belong2domain == 2) && (normal2Omega1_original(1)== -1))
            other_elem = elem_num + 1;
        elseif ((belong2domain == 2) && (normal2Omega1_original(1)== 1))
            other_elem = elem_num - 1; 
        % element in Omega1
        elseif ((belong2domain == 1) && (normal2Omega1_original(1)== -1))
            other_elem = elem_num - 1;
        elseif ((belong2domain == 1) && (normal2Omega1_original(1)== 1))
            other_elem = elem_num + 1; 
        else 
            error('the element normal do not match with the edge cut elements classification')
        end 
    else
        % arista coincides with vertical element boundary
        % element in Omega2
        if ((belong2domain == 2) && (normal2Omega1_original(2)== -1))
            other_elem = elem_num + nel_x;
        elseif ((belong2domain == 2) && (normal2Omega1_original(2)== 1))
            other_elem = elem_num - nel_x; 
        % element in Omega1
        elseif ((belong2domain == 1) && (normal2Omega1_original(2)== -1))
            other_elem = elem_num - nel_x;
        elseif ((belong2domain == 1) && (normal2Omega1_original(2)== 1))
            other_elem = elem_num + nel_x; 
        else 
            error('the element normal do not match with the edge cut elements classification')
        end 
    end 
    
    if ismember(other_elem,list1)
        list1_new = list1;
        list1_new(list1_new==other_elem) = [];
    else
        list1_new = list1; 
    end 

    if ismember(other_elem,list2)
        list2_new = list2;
        list2_new(list2_new==other_elem) = [];
    else
        list2_new = list2; 
    end 
    
    if isempty(list_cut) ~= 1
        if ismember(other_elem,list_cut(:,1))
            list_cut_new = list_cut; 
            list_cut_new(list_cut_new(:,1)==other_elem,:) = [];
        else
            list_cut_new = list_cut; 
        end 
    else
        list_cut_new = list_cut; 
    end 

    %warning('erase the elements from list_cut, list1, list2')

    % determine the intersections and normals
    Te = T(other_elem,1:4);
    LSe = InfoMesh.LS_new(Te); 
    % to obtain the coordinates and normals of the new element with a
    % coinciding interface at the arista, modify the level set values
    LSe(ismember(Te,T(elem_num,:))) = 0;
    Xe_corners = X(Te,:);
    [intersection1,intersection2] = findElemInter_intersection(Xe_corners,LSe);   % calculates the points (isop. coord.) of the intersection element-interface
    % normal pointing out of Omega1 considering the intersection points (interface linearly approx. in the element)
    normal2Omega1 = compute_normal_Omega1(intersection1,intersection2,Xe_corners,LSe);
end 
function list1_new = check_duplicated_classification(list1,list2)
    % check if there exist an element in list1 and list2 simultaneously
    a1 = ismember(list1(:,1),list2(:,1));
    list1_new = list1;
    list1_new(a1,:) = [];
end
function plot_crossed_elements(list_cut_elem,list_edge_elem1,list_edge_elem2,LABpoints,InfoMesh,InfoProblem)
    figure(120), clf
    X = InfoMesh.X; 
    T = InfoMesh.T;

    % Extract unique x and y values
    x_values = unique(X(:,1));  % X-values (horizontal)
    y_values = unique(X(:,2));  % Y-values (vertical)
    
    % Define the number of unique x points
    num_xpoints = length(x_values);
    num_ypoints = length(y_values);
    
    % Generate the grid where x varies along rows and y is constant for each row
    [x_grid1, y_grid1] = meshgrid(x_values, y_values);  % This creates the proper grid
    
    % Reshape the data in LS_new to fit this grid:
    % The number of rows corresponds to the length of y_values, and the number of columns
    % corresponds to the number of x_values.
    data2plot = reshape(InfoMesh.LS_new, [num_xpoints, num_ypoints])';
    
    % Now, plot the contour
    contourf(x_grid1, y_grid1, data2plot, 'FaceAlpha', 0.5);

    hold on
    
    % domain
    x1 = min(X(:,1));
    y1 = min(X(:,2));
    %z1 = min(X(:,3));
    x2 = max(X(:,1));
    y2 = max(X(:,2));
    %z2 = max(X(:,3));
    
    Xe = [x1 y1; x2 y1; x2 y2; x1 y2];
    addElem2(Xe,{'-k','linewidth',1})
    
    
    xplotLAB = LABpoints(1,:); %/(InfoProblem.L_ref/1000);
    yplotLAB = InfoProblem.L_ref*(InfoMesh.fin_y - InfoMesh.ini_y) - LABpoints(2,:);%(InfoProblem.L_ref/1000 + InfoLAB.LABy)/(InfoProblem.L_ref/1000);
    scatter(xplotLAB,yplotLAB,3,'Color',[0.8500 0.3250 0.0980],'Marker','s','MarkerEdgeColor',[0.8500 0.3250 0.0980]); %,'MarkerSize',1)
    
    plot([-1 0],[-1 -1],'b-')
    plot([-1 0],[-1 -1],'r-')
    plot([-1 0],[-1 -1],'g-')
    
    n_elem_tot = size(list_cut_elem,1);
    
    % elements interface correctly crossed
    for iCol = 1:n_elem_tot
      ce = list_cut_elem(iCol,1);
      for kEle = 1:length(ce)
         Te = T(ce(kEle),:);
         Xe = X(Te,:);
         addElem(Xe,{'-b'});
      end
    end
    
    n_ill_cross1 = size(list_edge_elem1,1);
    % elements interface correctly crossed
    for iCol = 1:n_ill_cross1
      ce = list_edge_elem1(iCol,1);
      for kEle = 1:length(ce)
         Te = T(ce(kEle),:);
         Xe = X(Te,:);
         addElem(Xe,{'-r'});
      end
    end
    
    n_ill_cross2 = size(list_edge_elem2,1);
    % elements interface correctly crossed
    for iCol = 1:n_ill_cross2
      ce = list_edge_elem2(iCol,1);
      for kEle = 1:length(ce)
         Te = T(ce(kEle),:);
         Xe = X(Te,:);
         addElem(Xe,{'-g'});
      end
    end
    
    legend('','','','elements correctly crossed','elements ill cross in \Omega_1','elements ill cross in \Omega_2','FontSize',12,'Location','SouthWest','AutoUpdate','off')
    xlabel('X axis','FontSize',12)
    ylabel('Y axis','FontSize',12)
    axis equal
    axis([x1 x2 y1 y2]);
end     