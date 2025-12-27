function plot_normal2interphase(vector,InfoMesh)
    % This function plots the normal vectors to the interphase lines
    % within a mesh. It uses the provided InfoMesh structure to extract
    % element and node information.
    % INPUT:
    %   vector    : Nx2 matrix where each row contains the normal vector
    %               components (normal_x, normal_y) for each interphase line.
    %   InfoMesh  : Structure containing mesh information including:
    %               - list_cut: Nx5 matrix with interphase line data
    %               - X: Node coordinates
    %               - T: Element connectivity matrix    
    % OUTPUT:
    %   A plot displaying the interphase lines and their normal vectors.

    % Sample data
    data = InfoMesh.list_cut; 
    
    % Assuming InfoMesh.X and InfoMesh.T are pre-loaded in your workspace
    % InfoMesh.X contains the coordinates of the nodes
    % InfoMesh.T contains the connectivity matrix (element to node relations)
    
    figure(120); %clf;
    hold on;
    axis equal;
    
    for i = 1:size(data, 1)
        % Extract data
        elem = data(i, 1); % Element number
        xi1 = data(i, 2);  % First isoparametric xi coordinate
        eta1 = data(i, 3); % First isoparametric eta coordinate
        xi2 = data(i, 4);  % Second isoparametric xi coordinate
        eta2 = data(i, 5); % Second isoparametric eta coordinate
        normal_x = vector(i,1);%data(i, 6); % Normal X-component
        normal_y = vector(i,2); % Normal Y-component
        
        % Extract element nodes from InfoMesh.T
        element_nodes = InfoMesh.T(elem, :); % Connectivity for the current element
        Xe = InfoMesh.X(element_nodes, :); % Coordinates of the element's nodes
        
        % Find the physical coordinates of the intersection points
        points1 = [xi1 eta1; xi2 eta2];
        [N_point,~,~] = shapeFunctions(InfoMesh.elemType,InfoMesh.nne,points1);
        point1 = N_point(1,:)*Xe;
        point2 = N_point(2,:)*Xe;
        
        % Plot the interface line
        plot([point1(1), point2(1)], [point1(2), point2(2)], 'b-', 'LineWidth', 2);
        
        % Midpoint of the interface line
        midpoint = (point1 + point2) / 2;
        
        % Plot the normal vector
        quiver(midpoint(1), midpoint(2), normal_x, normal_y, 8000, 'r', 'LineWidth', 2);
    end
    
    title('Element Interfaces and Normals');
    xlabel('X');
    ylabel('Y');
    grid on;
    hold off;
end 