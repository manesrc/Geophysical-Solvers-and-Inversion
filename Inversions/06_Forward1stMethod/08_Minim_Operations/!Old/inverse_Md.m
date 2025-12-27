function [invMd] = inverse_Md(UnscaledM,selected_method,tolerance)
% function to compute the non-unique inverse matrix relating the position
% of the interface to the degrees of freedom of the crossed elements
% INPUTS
% UnscaledM: mass matrix computed for Nitsche's method, M = ∫ N_i * N_j d\Gamma [Unscaled because it doesn't consider ß]
% selected_method: two cases,  1) Lumped mass [sum the rows converts the matrix to a diagonal and inverts]
%                                                   CODE: selected_method ==  0
%                                               2) Pseudo_inverse: M = U*S*V', inverse is M^(-1) = V*S^(-1)*U', and the S^(-1) is considered up to a certain tolerance for the eigenvalues, tolerance
%                                                   CODE: selected_method ==  1
% tolerance: minimum eigenvalue to be consired for the Pseudo_inverse case


if selected_method == 0
    invMd = zeros(size(UnscaledM));
    for pp = 1:size(UnscaledM,1)
        invMd(pp,pp) = 1/sum(UnscaledM(pp,:));
    end 
else
    [U,S,V] = svd(full(UnscaledM)); % svd of UnscaledM, then apply the tolerance
    S_tol_inv = zeros(size(S));       % new matrix to store eigenvalues bigger than tol.
    for ii =1:size(S,1)
        if abs(S(ii,ii)) < tolerance
            S_tol_inv(ii,ii) = 0;           % replace eig 
        else
            S_tol_inv(ii,ii) = 1/S(ii,ii);    % keep eig
        end 
    end 
    invMd = V*S_tol_inv*U'; 
end 

end 