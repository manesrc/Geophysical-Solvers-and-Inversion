function U_b_sym = generate_symmetric_base(U_basis,nel_x,nel_y)
% generate_symmetric_base - Reduces the original basis to its half by summing the
% symmetric components with respect to the Y-Y axis.
% INPUTS:
%   U_basis - Matrix of basis velocity components (size: m x n)
%   nel_x - Number of elements in the X direction (must be even)
%   nel_y - Number of elements in the Y direction
%
% OUTPUT:
%   U_b_sym - Symmetric basis matrix

    % Validate that nel_x is even
    if mod(nel_x, 2) ~= 0
        error('nel_x must be an even number for symmetric reduction.');
    end
    U_b_sym = zeros(size(U_basis,1), size(U_basis,2)/2);
    % Remember U_basis is constructed from top to bottom in Y-Direction

    for ii = 1:nel_y
        ini_left = (ii-1)*nel_x+1;
        fin_left = ini_left + nel_x/2 - 1;
        left_side = ini_left:1:fin_left;
        ini_right = fin_left+1;
        fin_right = ii*nel_x; 
        right_side = flip(ini_right:1:fin_right);
        U_new = U_basis(:,left_side) + U_basis(:,right_side);
        ind2save = (ii-1)*(nel_x/2)+1:ii*(nel_x/2);
        U_b_sym(:,ind2save) = U_new; 
    end 
end