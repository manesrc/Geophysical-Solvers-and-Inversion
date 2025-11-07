function [sortedX,sortedY] = sortLABPoints2plot(X,Y,fig_num,plot2)
% This function sorts the LAB interface points in a parametric manner for plotting
% OUTPUT:
%   sortedX, sortedY: Parametrically sorted coordinates of the LAB interface points
% INPUT:
%   X, Y: Coordinates of the LAB interface points
%   fig_num: Figure number for plotting
%   plot2: Flag to indicate whether to plot the sorted points (1 = plot,
%          0 = do not plot)
    % check point uniqueness
    points = [X(:), Y(:)];

    % Remove duplicate rows (points)
    [uniquePoints, uniqueIdx] = unique(points, 'rows', 'stable');
    
    % Extract the unique X and Y points
    uniqueX = uniquePoints(:, 1);
    uniqueY = uniquePoints(:, 2);

    % Number of points
    nPoints = length(uniqueX);
    
    % Create an array to hold the ordered points
    sortedX = zeros(nPoints, 1);
    sortedY = zeros(nPoints, 1);
    
    % Track which points have been visited
    visited = false(nPoints, 1);
    
    % Start at the first point
    [~,idx1] = sort(uniqueX,'ascend');
    currentIdx = idx1(1);
    sortedX(1) = uniqueX(currentIdx);
    sortedY(1) = uniqueY(currentIdx);
    visited(currentIdx) = true;
    
    % Iterate to find the next nearest point and build the path
    for i = 2:nPoints
        minDist = inf;
        nextIdx = -1;
        
        % Find the nearest unvisited point
        for j = 1:nPoints
            if ~visited(j)
                dist = sqrt((uniqueX(currentIdx) - uniqueX(j))^2 + (uniqueY(currentIdx) - uniqueY(j))^2);
                if dist < minDist
                    minDist = dist;
                    nextIdx = j;
                end
            end
        end
        
        % Update the sorted points list and mark the point as visited
        currentIdx = nextIdx;
        sortedX(i) = uniqueX(currentIdx);
        sortedY(i) = uniqueY(currentIdx);
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