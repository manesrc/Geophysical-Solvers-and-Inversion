function U_basis_ordered = reorder_Ubasis(Ubasis,nel_x,nel_y)
    % this function aims at re-ordering the originally computed Ubasis (that goes 
    % from top-left to bottom right) as the mesh numbering that goes from bottom 
    % left to top right. 
    % OUTPUT: 
    % U_basis_ordered: reordered Ubasis
    % INPUT:
    % Ubasis: original velocity basis
    % nel_x: number of elements in x direction
    % nel_y: number of elements in y direction
    
    U_basis_ordered = zeros(size(Ubasis)); 
    num_of_comp = size(Ubasis,2);
    % ini_original_basis = nel_x*nel_y - inv_nelx * nel_y + 1
    % fin_original_basis = nel_x*nel_y - (inv_nelx-1) * nel_y
    % inv_nelx = nel_x:-1:1;
    inv_nelx = nel_x:-1:1;
    for ii = 1:nel_x
        inv_nelx_ii = inv_nelx(ii);
        comp_col_ii = flip(ii:nel_x:num_of_comp);
        for jj = 1:nel_y
            ini_original_basis = (nel_x - inv_nelx_ii)*nel_y + jj;
            comp_velo = Ubasis(:,ini_original_basis);
            U_basis_ordered(:,comp_col_ii(jj)) = comp_velo;
        end 
    end
end 