function [list_Omega1,list_Omega2,list_cut_elem,list_edge_elem1,list_edge_elem2] = findElemNitsche_LS(InfoMesh,InfoLAB,tol,plot1)
% Level set
LS = InfoMesh.LS_mesh;
% mesh matrices
X = InfoMesh.X;
T = InfoMesh.T;

% data to save
list_Omega1 = [];
list_Omega2 = [];
list_cut_elem = [];
list_edge_elem1 = [];
list_edge_elem2 = [];


for jj = 1:size(T,1)
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

    if count_neg == count_tot || (count_neg == count_tot-1 && count_zeros == 1)
        list_Omega2 = [list_Omega2; jj];      % element in Omega 2
    elseif count_pos == count_tot || (count_pos == count_tot-1 && count_zeros == 1)
        list_Omega1 = [list_Omega1; jj];      % element in Omega 1
    elseif count_pos > 0 && count_neg > 0
        % intersection points
        Xe_corners = X(Te,:);   % X,Y points of the corners of the element
        [intersec_1,intersec_2] = findElemInter_intersection(Xe_corners,LSe);   % calculates the points (isop. coord.) of the intersection element-interface
        % normal pointing out of Omega1 considering the intersection points (interface linearly approx. in the element)
        normal2Omega1 = compute_normal_Omega1(intersec_1,intersec_2,Xe_corners,LSe);
        % list of element crossed by interface + geometric characteristics of interest
        list_cut_elem = [list_cut_elem; jj intersec_1 intersec_2 normal2Omega1];      
    elseif (count_zeros == 2) && (count_pos == 2 || count_neg == 2)      % case line coinciding with arista
        % intersection points (coincident with zeros)
        Xe_corners = X(Te,:);   % X,Y points of the corners of the element
        [intersec_1,intersec_2] = findElemInter_intersection(Xe_corners,LSe);   % calculates the points (isop. coord.) of the intersection element-interface
        % normal pointing out of Omega1 considering the intersection points (interface linearly approx. in the element)
        normal2Omega1 = compute_normal_Omega1(intersec_1,intersec_2,Xe_corners,LSe);
        if count_neg == 2
            list_edge_elem2 = [list_edge_elem2; jj intersec_1 intersec_2 normal2Omega1];
        else
            list_edge_elem1 = [list_edge_elem1; jj intersec_1 intersec_2 normal2Omega1];
        end 
    else
        error('Si entra acá: FALTAN COSAS')
    end 
end 

%% plot crossed elements 

if plot1 == 1
    figure(120), clf
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
    
    
    xplotLAB = InfoLAB.LABx/(InfoLAB.maxDepth/1000);
    yplotLAB = (InfoLAB.maxDepth/1000 + InfoLAB.LABy)/(InfoLAB.maxDepth/1000);
    plot(xplotLAB,yplotLAB,'Color',[0.8500 0.3250 0.0980],'Marker','s','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize',1)
    
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
    
    legend('','','elements correctly crossed','elements ill cross in \Omega_1','elements ill cross in \Omega_2','FontSize',12,'Location','SouthWest','AutoUpdate','off')
    
    xlabel('X axis','FontSize',12)
    ylabel('Y axis','FontSize',12)
    axis equal
    axis([0 1 0 1]);
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

%% plot crossed elements 
%figure(120), clf
% plot over Level-Set
% figure(1)
% hold on
% 
% % domain
% x1 = min(X(:,1));
% y1 = min(X(:,2));
% %z1 = min(X(:,3));
% x2 = max(X(:,1));
% y2 = max(X(:,2));
% %z2 = max(X(:,3));
% 
% Xe = [x1 y1; x2 y1; x2 y2; x1 y2];
% addElem2(Xe,{'-k','linewidth',1})
% 
% 
% xplotLAB = InfoLAB.LABx/(InfoLAB.maxDepth/1000);
% yplotLAB = (InfoLAB.maxDepth/1000 + InfoLAB.LABy)/(InfoLAB.maxDepth/1000);
% plot(xplotLAB,yplotLAB,'Color',[0.8500 0.3250 0.0980],'Marker','s','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerSize',1)
% 
% plot([-1 0],[-1 -1],'b-')
% plot([-1 0],[-1 -1],'r-')
% plot([-1 0],[-1 -1],'g-')
% 
% n_elem_tot = size(elements_interface,1);
% 
% % elements interface correctly crossed
% for iCol = 1:n_elem_tot
%   ce = elements_interface(iCol,1);
%   for kEle = 1:length(ce)
%      Te = T(ce(kEle),:);
%      Xe = X(Te,:);
%      addElem(Xe,{'-b'});
%   end
% end
% 
% n_ill_cross1 = size(elem_illCross_1,1);
% % elements interface correctly crossed
% for iCol = 1:n_ill_cross1
%   ce = elem_illCross_1(iCol,1);
%   for kEle = 1:length(ce)
%      Te = T(ce(kEle),:);
%      Xe = X(Te,:);
%      addElem(Xe,{'-r'});
%   end
% end
% 
% n_ill_cross2 = size(elem_illCross_2,1);
% % elements interface correctly crossed
% for iCol = 1:n_ill_cross2
%   ce = elem_illCross_2(iCol,1);
%   for kEle = 1:length(ce)
%      Te = T(ce(kEle),:);
%      Xe = X(Te,:);
%      addElem(Xe,{'-g'});
%   end
% end
% 
% legend('','','','','elements correctly crossed','elements ill-crossed in \Omega_1','elements ill-crossed in \Omega_2','Location','SouthWest')
% 
% xlabel('X axis')
% ylabel('Y axis')
% axis equal
% axis([0 1 0 1]);