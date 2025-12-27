function [U_squares, components] = sum_square_velocities(U_basis,nel_x,nel_y)
% sum_square_velocities - Reduces the original basis to its quarter by
% summing two elements in X- and two elements in Y-directions
% INPUTS:
%   U_basis - Matrix of basis velocity components (size: m x n)
%   nel_x - Number of elements in the X direction (must be even)
%   nel_y - Number of elements in the Y direction
% OUTPUT:
%   U_squares - Symmetric basis matrix


    % Validate that nel_x and nel_y are even
    if mod(nel_x, 2) ~= 0 || mod(nel_y, 2) ~= 0
        error('nel_x and nel_y must be even numbers for square reduction.');
    end
    
    % Initialize the reduced basis matrix
    U_squares = zeros(size(U_basis, 1), (size(U_basis, 2) / 4));
    components = zeros((size(U_basis, 2) / 4), 4);

    % Counter to keep track of the column index in the reduced matrix
    counter = 0;

    % Loop over half the elements in the Y direction
    for jj = 1:(nel_y / 2)
        % Loop over half the elements in the X direction
        for ii = 1:(nel_x / 2)
            % Calculate indices for the four components in the original matrix
            base_index = (jj - 1) * 2 * nel_x + (ii - 1) * 2 + 1;
            base_ind1 = base_index;
            base_ind2 = base_ind1 + 1;
            vert_ind1 = base_ind1 + nel_x;
            vert_ind2 = base_ind2 + nel_x;
          %  [base_ind1 base_ind2 vert_ind1 vert_ind2]

            % Sum the four components
            U_new = U_basis(:, vert_ind1) + U_basis(:, vert_ind2) + ...
                    U_basis(:, base_ind1) + U_basis(:, base_ind2);

            % Increment the counter and save the new summed component
            counter = counter + 1;
            U_squares(:, counter) = U_new;
            components(counter, : ) = [base_ind1 base_ind2 vert_ind1 vert_ind2]; 
        end
    end
end
