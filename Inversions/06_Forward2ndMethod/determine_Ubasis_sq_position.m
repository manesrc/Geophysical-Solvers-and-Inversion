function determine_Ubasis_sq_position(nel_x,nel_y)
% this function determines which elements in the original basis goes to
% which components of the basis grouped by four
    num_comp_tot = nel_x*nel_y; 
    num_comp_Sq_basis = num_comp_tot/4; 
    comp_in_sq_basis = zeros(num_comp_Sq_basis,4);
    
    for ii = 1:num_comp_Sq_basis
        comp_in_sq_basis(ii,:) = [ii ii+1 ]
    end 


end 