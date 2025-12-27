function [sortedX,sortedY] = sortLABPoints2plot(X,Y,fig_num,plot2)
% sortLABPoints2plot - Sorts the LAB 2D points (X,Y) in a parametric order for plotting.
% Inputs:
%    X      - Array of X coordinates of LAB points
%    Y      - Array of Y coordinates of LAB points
%    fig_num - Figure number for plotting
%    plot2   - Flag to indicate whether to plot the sorted points (1 = plot, 0 = do not plot)
% Outputs:
%    sortedX - Array of sorted X coordinates
%    sortedY - Array of sorted Y coordinates

    % Number of points
    nPoints = length(X);
    
    % Create an array to hold the ordered points
    sortedX = zeros(nPoints, 1);
    sortedY = zeros(nPoints, 1);
    
    % Track which points have been visited
    visited = false(nPoints, 1);
    
    % Start at the first point
    [~,idx1] = sort(X,'ascend');
    currentIdx = idx1(1);
    sortedX(1) = X(currentIdx);
    sortedY(1) = Y(currentIdx);
    visited(currentIdx) = true;
    
    % Iterate to find the next nearest point and build the path
    for i = 2:nPoints
        minDist = inf;
        nextIdx = -1;
        
        % Find the nearest unvisited point
        for j = 1:nPoints
            if ~visited(j)
                dist = sqrt((X(currentIdx) - X(j))^2 + (Y(currentIdx) - Y(j))^2);
                if dist < minDist
                    minDist = dist;
                    nextIdx = j;
                end
            end
        end
        
        % Update the sorted points list and mark the point as visited
        currentIdx = nextIdx;
        sortedX(i) = X(currentIdx);
        sortedY(i) = Y(currentIdx);
        visited(currentIdx) = true;
    end
    
    % Plot the sorted points to verify the result
    if plot2 ==1
        figure(fig_num);
        plot(sortedX, sortedY, '-o');
        title('Parametrically Sorted Points');
        xlabel('X');
        ylabel('Y');
        grid on;
    end 
end 