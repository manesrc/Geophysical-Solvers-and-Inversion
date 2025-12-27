function [prior_bounds,relev_DOF] = generate_prior(LABmin, LABmax, ZonesLABPrior, X, maxDepth, plot_idea)
    
    % Generate bounds for multiple areas
    [handling_functions,releDOF] = compute_handling_functions(LABmin, LABmax, ZonesLABPrior, X, maxDepth, plot_idea);

    % save the modified DOF to a vector
    
    % Ensure all DOFs are column vectors
    releDOF = cellfun(@(x) x(:), releDOF, 'UniformOutput', false);
    % Concatenate them into a single vector
    relev_DOF = unique(cell2mat(releDOF));
    
    % all the DOFs have priors
    prior_bounds = zeros(size(X,1),2);
    prior_bounds(:,1) = -Inf;
    prior_bounds(:,2) = Inf;

    % modify those in the specific areas
    for i = 1:length(handling_functions)
        % Get segment functions and range
        LB = handling_functions(i).LB;
        UB = handling_functions(i).UB;
        
        % Compute bounds for the relevant nodes
        DOFinRelevArea = releDOF{i};    % get the DOFs to condition
        posX = X(DOFinRelevArea,1);     % compute X-coord
        posY = X(DOFinRelevArea,2);     % compute Y-coord
        % compute the values of the lower and upper prior functions
        segment_bounds = [UB(posX, posY), LB(posX, posY)]; 
        % save
        prior_bounds(DOFinRelevArea,:) = segment_bounds; 

        % exhaustive checks!!
        % figure(1); hold on; 
        % scatter3(posX,posY,0*posX,'filled'); 
        % scatter3(posX,posY,LB(posX,posY),'filled'); 
        % scatter3(posX,posY,UB(posX,posY),'filled'); 
    end

end

function [handling_functions,defRelevDOFperFunc] = compute_handling_functions(vectMinimumYDepth, vectMaximumYDepth, vectBoundariesX, X, maxDepth, plot_idea)
    % Input: 
    % vectMinimumY: Vector of minimum Y values for each area
    % vectMaximumY: Vector of maximum Y values for each area
    % vectBoundariesX: Boundaries in X defining transitions between areas
    % X: Mesh points matrix
    % maxDepth: Maximum depth for the domain
    % plot_idea: Boolean to enable/disable plotting

    % check inputs make sense
    if any(diff([min(X(:,1)), vectBoundariesX, max(X(:,1))]) <= 0)
        error('Overlap or gaps detected in boundary definitions.');
    end

    % Number of segments
    nSegments = length(vectMinimumYDepth);

    % Add start (0) and end (max(X)) boundaries
    boundaries = [min(X(:,1)), vectBoundariesX, max(X(:,1))];
    
    % Preallocate function handles
    handling_functions = struct('LB', {}, 'UB', {}, 'rangeX', {});

    % save DOF relevant for each area of the domain
    all_DOF = 1:1:size(X,1); 
    defRelevDOFperFunc = cell(nSegments,1);


    % Iterate over each segment
    for i = 1:nSegments
        LABminDepth = vectMinimumYDepth(i);
        LABmaxDepth = vectMaximumYDepth(i);
        
        % Segment-specific bounds
        phi_UB_at_00 = -(maxDepth - LABminDepth);    % values of phi for the UB function in (x,y) = (0,0)
        phi_LB_at_00 = -(maxDepth - LABmaxDepth);    % values of phi for the LB function in (x,y) = (0,0)
        
        % Linear functions for UB and LB
        mx = 0; my = 1; % Assuming vertical gradients
        LB = @(x, y) phi_LB_at_00 + mx*x + my*y;
        UB = @(x, y) phi_UB_at_00 + mx*x + my*y;

        % Store functions and range
        handling_functions(i).LB = LB;
        handling_functions(i).UB = UB;
        handling_functions(i).rangeX = [boundaries(i), boundaries(i+1)];

        % Define DOFs relevant for each domain area
        defRelevDOFperFunc{i}= all_DOF(X(:,1) >= boundaries(i) & X(:,1) <= boundaries(i+1) & ...
            X(:,2) >= -phi_LB_at_00 & X(:,2) <= -phi_UB_at_00);
    end

    % Optional: Plot the segments
    if plot_idea
        figure; hold on;
        for i = 1:nSegments
            % Create mesh for plotting
            x_range = linspace(handling_functions(i).rangeX(1), handling_functions(i).rangeX(2), 20);
            y_range = linspace(min(X(:,2)), max(X(:,2)), 20);
            [x2plot, y2plot] = meshgrid(x_range, y_range);

            % Plot upper and lower bounds
            surf(x2plot, y2plot, handling_functions(i).UB(x2plot, y2plot), 'FaceColor', 'red', 'FaceAlpha', 0.5);
            surf(x2plot, y2plot, handling_functions(i).LB(x2plot, y2plot), 'FaceColor', 'blue', 'FaceAlpha', 0.5);
        end
        xlinspac1 = linspace(min(X(:,1)), max(X(:,1)), 20);
        ylinspac1 = linspace(min(X(:,2)), max(X(:,2)), 20);
        [x2plot1,y2plot1] = meshgrid(xlinspac1,ylinspac1);
        surf(x2plot1,y2plot1,0*x2plot1)
        xlabel('X-Axis [m]'); ylabel('Y-Axis [m]'); zlabel('\phi-value [m]');
        legend('UB', 'LB','','','','','plane', 'Location', 'Best');
        title('Upper and Lower Bounds for Each Segment');
        view(3)
    end
end