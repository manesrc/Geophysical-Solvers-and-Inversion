function [U_basis_new, flag1] = update_basis_data(U_basis, ...
    components_velo_basis, Omega1_new, Omega1_old)
% Update the velocity basis based on changes in Omega1
% INPUTS:
% - U_basis: Current velocity basis matrix
% - components_velo_basis: Groups of velocity components in the basis
% - Omega1_new, Omega1_old: Sets of elements in the lithosphere for current and previous iteration
% OUTPUTS:
% - U_basis_new: Updated velocity basis matrix
% - flag1: Logical array indicating changed components (1 if updated, 0 otherwise)

    % Determine components to trim for the new and old Omega1
    components2trim_new = logical(sum(ismember(components_velo_basis, Omega1_new), 2) == 4);
    if isempty(Omega1_old)
        components2trim_old = false(size(components2trim_new));
    else
        components2trim_old = logical(sum(ismember(components_velo_basis, Omega1_old), 2) == 4);
    end

    % Update U_basis by removing trimmed components
    U_basis_new = U_basis;
    U_basis_new(:, components2trim_new) = []; 

    % Compute the flag indicating changed components
    flag1 = [components2trim_new components2trim_old];  % 1 for updated, 0 otherwise
end
