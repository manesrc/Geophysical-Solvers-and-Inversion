function plotBothSub_dim_EGU(temp1,temp2,InfoMesh,InfoLAB,InfoProblem,nodes1,nodes2,grad2plot,xplotgrad,fig_num,varargin)
% unfold mesh
T = InfoMesh.T;
X = InfoMesh.X;

ne_x = InfoMesh.nel_x; % number of elements in X (horizontal) direction
ne_y = InfoMesh.nel_y; % number of elements in Y (vertical) direction

%T_LAB = InfoProblem.T_LAB * InfoProblem.T_ref; % [K]
T_0_abs = 273; % [K]

temp_inC1 = temp1 - T_0_abs; % Temp in K to ºC
temp_inC2 = temp2 - T_0_abs; % Temp in K to ºC



%% interpolate the solution
% define where in each element will interpolate the solution
% if ne_x >= 100
%     isop = linspace(-0.999,0.999,10)';
% else
%     isop = linspace(-0.999,0.999,50)';
% end 
isop = linspace(-0.999,0.999,10*(400/ne_x))';
%isop = linspace(-0.999,0.999,50)';
[chi_isop,eta_isop] = meshgrid(isop,isop);
chi_isop2 = reshape(chi_isop,size(chi_isop,1)*size(chi_isop,2),1);
eta_isop2 = reshape(eta_isop,size(eta_isop,1)*size(eta_isop,2),1);

% generate to fake sub-domains
% base elements
be = (1:1:ne_x)';
% top elements
te = ((ne_x*ne_y-ne_x)+1:1:ne_x*ne_y)';

% to determine how many elements to consider in each plot, find maximums
% and minimums of the interface
minLAB = min(InfoLAB.maxDepth*(InfoMesh.fin_y - InfoMesh.ini_y) + InfoLAB.LABy*1000)/(InfoLAB.maxDepth);
maxLAB = max(InfoLAB.maxDepth*(InfoMesh.fin_y - InfoMesh.ini_y) + InfoLAB.LABy*1000)/(InfoLAB.maxDepth);
elem2add_Omega1 = ceil(minLAB/( (InfoMesh.fin_y-InfoMesh.ini_y) /ne_y));
elem2add_Omega2 = ceil(maxLAB/( (InfoMesh.fin_y-InfoMesh.ini_y) /ne_y));
elem_perCol_Omega1 = length(elem2add_Omega1:ne_y);
fake_Omega1 = zeros(ne_x*elem_perCol_Omega1,1);
elem_perCol_Omega2 = length(1:elem2add_Omega2);
fake_Omega2 = zeros(ne_x*elem_perCol_Omega2,1);
for ii_fake = 1:length(be)
    columnElements = be(ii_fake,1):ne_x:te(ii_fake,1);
    addOmega1 = columnElements(elem2add_Omega1:end);
    count_ini1 = (ii_fake-1)*elem_perCol_Omega1+1; count_fin1 = ii_fake*elem_perCol_Omega1;
    fake_Omega1(count_ini1:count_fin1,1) = addOmega1';
    addOmega2 = columnElements(1:elem2add_Omega2);
    count_ini2 = (ii_fake-1)*elem_perCol_Omega2+1; count_fin2 = ii_fake*elem_perCol_Omega2;
    fake_Omega2(count_ini2:count_fin2,1) = addOmega2';
end 

%% select the subdomain
problem_int1 = [1,2];
for ii_plot = 1:2
    problem_int = problem_int1(ii_plot);
    if problem_int == 1      
        if isempty(InfoMesh.list_edge1) ~= 1 && isempty(InfoMesh.list_cut) ~=1
            omeguita1 = [InfoMesh.list1; InfoMesh.list_cut(:,1); InfoMesh.list_edge1(:,1)];
            elements_interface = [InfoMesh.list_cut(:,1); InfoMesh.list_edge1(:,1)];
        elseif isempty(InfoMesh.list_edge1) == 1
            omeguita1 = [InfoMesh.list1; InfoMesh.list_cut(:,1)];
            elements_interface = [InfoMesh.list_cut(:,1)];
        elseif isempty(InfoMesh.list_cut) == 1
            omeguita1 = [InfoMesh.list1; InfoMesh.list_edge1(:,1)];
            elements_interface = [InfoMesh.list_edge1(:,1)];
        else
            error('Problemas')
        end 
        num_of_elements = size(fake_Omega1,1);
        domain_interest = fake_Omega1;
        omeguita_interest = omeguita1;
        kk_el = elem_perCol_Omega1;
        X_isop = [eta_isop2 chi_isop2];
        [N,~,~] = shapeFunctions(InfoMesh.elemType,InfoMesh.nne,X_isop);
        temp = temp_inC1;
        nodes_int = nodes1;
    else
        if isempty(InfoMesh.list_edge2)~= 1 && isempty(InfoMesh.list_cut) ~= 1
            omeguita2 = [InfoMesh.list2; InfoMesh.list_cut(:,1); InfoMesh.list_edge2(:,1)];
            elements_interface = [InfoMesh.list_cut(:,1); InfoMesh.list_edge2(:,1)];
        elseif isempty(InfoMesh.list_edge2) == 1
            omeguita2 = [InfoMesh.list2; InfoMesh.list_cut(:,1)];
            elements_interface = [InfoMesh.list_cut(:,1)];
        elseif isempty(InfoMesh.list_cut) == 1
            omeguita2 = [InfoMesh.list2; InfoMesh.list_edge2(:,1)];
            elements_interface = [InfoMesh.list_edge2(:,1)];
        else
            error('Problemas')
        end 
        nel_per_column = elem2add_Omega2;
        num_of_elements = size(fake_Omega2,1);
        domain_interest = fake_Omega2;
        omeguita_interest = omeguita2;
        kk_el = min(nel_per_column);
        X_isop1 = [eta_isop2 chi_isop2];
        X_isop = X_isop1;
        X_isop(:,2) = X_isop1(:,2)*(-1);
        [N,~,~] = shapeFunctions(InfoMesh.elemType,InfoMesh.nne,X_isop);
        temp = temp_inC2;
        nodes_int = nodes2;
    end 
    
    nNodes2plot_iso = size(X_isop,1);
    Xe2plot = zeros(num_of_elements*nNodes2plot_iso,2);
    
    % plot matrix T_r
    n_grid = size(isop,1); % the points in the grid are a quadrilateral (n_grid_x = n_grid_y)
    T_r = zeros(n_grid*kk_el,n_grid*ne_x);
    % tolerance
    tolerance = InfoProblem.tolerance;
    
    
    for ii = 1:num_of_elements
        elem_interest = domain_interest(ii);
        Te = T(elem_interest,:);
        Xe = X(Te,:);
        % temperatures in nodes to interpolate
        temp_nodes_element = zeros(size(T,2),1);
        LSe = zeros(length(Te),1);
        if ismember(elem_interest,omeguita_interest) % temperatures relevant to the plot
            for jj = 1:size(T,2)
                %nodal temperatures
                temp_nodes_element(jj) = temp(nodes_int == Te(jj));
                % find the elemental Level set values
                is_bigger_than_tolerance = abs(InfoMesh.LS_mesh(Te(jj)))>tolerance;
                LSe(jj) = (is_bigger_than_tolerance) * InfoMesh.LS_mesh(Te(jj));
            end 
            % check if the element is in the interface area
            if ismember(elem_interest,elements_interface)
                % Gauss weight with shape functions
                He = N*LSe;
                if problem_int == 1
                    He(He < -tolerance) = NaN;
                else
                    He(He > tolerance) = NaN;
                end 
                He(~isnan(He)) = 1;
            else
                He = ones(size(N,1),1);
            end 
        else
            temp_nodes_element = NaN*ones(size(T,2),1);
            He = zeros(size(N,1),1);
        end 
        % points X,Y of the temperature
        Xe2plot((1+(ii-1)*nNodes2plot_iso:ii*nNodes2plot_iso),:) = N*Xe;   
        % temperatures 
        Temps = (N*temp_nodes_element).*He;     
        % elemental matrix of temperature results
        T_r_el = reshape(Temps,[n_grid n_grid])';
        % fill the general  matrix of results
        iniCol = (ceil(ii/kk_el)-1)*n_grid+1; 
        endCol = ceil(ii/kk_el)*n_grid;

        if floor(ii/kk_el) ~= ii/kk_el
            move_index1 = (ii/kk_el - floor(ii/kk_el))*kk_el;
        else
            move_index1 = kk_el;
        end 
        if problem_int == 1
            iniRow = int16((move_index1-1)*n_grid + 1);
            endRow = int16(move_index1*n_grid);
        else
            endRow = int16((kk_el - (move_index1-1))*n_grid);
            iniRow = int16((kk_el - move_index1)*n_grid + 1);
        end 
            
        T_r(iniRow:endRow,iniCol:endCol) = T_r_el;
    end 
    
    %% Plot definition matrices
    xu = unique(round(Xe2plot(:,1),6));
    yu = 1 - unique(round(Xe2plot(:,2),6)); % write in depth: L_max - x & in km /1000

    % Ensure that xu and yu have the same dimensions as T_r
    xu = xu(1:size(T_r, 2)); % Match the number of columns
    yu = yu(1:size(T_r, 1)); % Match the number of rows

    fonti_size = 16;

    if ii_plot == 1
        xu1 = xu;
        yu1 = yu;
        T_r1 = T_r;
        figure(fig_num);
        clf
        hold on
        sf1 = subplot(2,1,1);
        sf2 = subplot(2,1,2);
        sf1.Position = [0.1300    0.250    0.7750    0.700];
        %sf2.Position = [0.1300    0.1100    0.7750    0.1100];
        sf2.Position = [0.1300    0.1100    0.7750    0.080];
        sf2.Box = 'on';
        sf1.Box = 'on';
        hold on
        sf1.XGrid = 'on';
        sf1.YGrid = 'on';

        %%
        hold(sf1,'on')
        %plot(sf1,Xinterf/1000,Yinterf,'b--','LineWidth',1.5)
        set(sf1,'FontSize',fonti_size)
        %sf1.XTick = [];
        sf1.XTickLabel = [];
        ylabel(sf1,'Depth [km]','FontSize',fonti_size)
        %grad2plot2= [grad2plot; grad2plot(end)];
        %plot(sf2,(InfoProblem.L_ref/1000)*xplotgrad,grad2plot_mean','k-','LineWidth',1.0)
        vect1 = grad2plot;
        plot(sf2,(InfoProblem.L_ref/1000)*xplotgrad,grad2plot','b-','LineWidth',1.0)
        axis_info = [0 (InfoProblem.L_ref/1000) min(min(vect1)) max(max(vect1))];
        axis(axis_info);% round(1.1*max(1000*grad2plot),1)]);
        hold(sf2,'on')
        sf2.XGrid = 'on';
        sf2.YGrid = 'on';
        sf2.XLimitMethod = 'tight';
        %sf2.YLim = [0.2 0.6];
        %sf2.YLimitMethod = 'tight';
        %yticks(sf2,[0.2 0.4 0.6]) 

        sf2.YLim = [-5 5];
        sf2.YTick = -5.0:5.0:5.0;
        sf2.YTickLabel = {'-5.0','0.0','5.0'};
        xlabel(sf2,'X [km]','FontSize',fonti_size)
        ylabel(sf2,'\nabla T [ºC/km]','FontSize',fonti_size)
        set(sf2,'FontSize',fonti_size)
    else
        yu = sort(yu,'ascend'); % for depth 'ascend' for "distance to the surf" put 'descend'
        figure(fig_num)
        hold(sf1,'on')
        % Set color limits
        %contours2 = [1100:100:1500 T_LAB 1600:20:1700 1700:100:2000]' - InfoProblem.T_sup * InfoProblem.T_ref;
        contours2 = [1000:100:1200 1300:20:1360 1400:100:1600]';
        [C2,h2] = contourf(sf1,(InfoProblem.L_ref/1000)*xu,(InfoProblem.L_ref/1000)*yu,T_r,contours2);     
        %contours1 = [293; 500; 800; 1100; 1300; 1500; T_LAB] - InfoProblem.T_sup * InfoProblem.T_ref;
        contours1 = [20 200:300:800 1100 1300 1320];
        [C1,h1] = contourf(sf1,(InfoProblem.L_ref/1000)*xu1,(InfoProblem.L_ref/1000)*yu1,T_r1,contours1);
        iso_line_width = 1.2;
        limits = [20 1600];
        colormap(jet)
        caxis(sf1, limits);
        % create a colorbar with custom tick levels 
        contour_mod = colorbar; %contourcbar('eastoutside');
        caxis([20 1600]);
        contour_mod.XLabel.String = 'Temperature [ºC]';
        %contour_mod.Ticks = levels1;
        contour_mod.Position = [sf1.Position(1)+0.7 sf1.Position(2) 0.05 0.7];
        contour_mod.FontSize = 14;
        %contour_mod.TickLength = 0.09;
        % reshape the figure 
        daspect(sf1,[1 1 1])
        sf1.YDir = 'reverse'; 
        %contour_mod = contourcmap('jet',cmap,'Colorbar','on','TitleString','Temperature [K]','FontSize',14);
        set(sf2,'Position',[0.255 0.11 0.525 0.11])
        %colormap('jet')
        h1.LineWidth = iso_line_width;
        h2.LineWidth = iso_line_width;
        %v2 = [1000 1200 1400 T_LAB:20:T_LAB+60 1700:100:1800];
        clabel(C1,h1,'FontSize',14);
        clabel(C2,h2,contours2,'FontSize',14);
        Xinterf = InfoLAB.LABx;%/(InfoLAB.maxDepth/1000); 
        %Yinterf = (InfoProblem.L_ref/1000)+InfoLAB.LABy;
        plot(sf1,Xinterf,abs(InfoLAB.LABy),'k-','LineWidth',iso_line_width);
        plot(sf1,Xinterf,abs(InfoLAB.LABy),'b--','LineWidth',3.0)
    end 



end 



end