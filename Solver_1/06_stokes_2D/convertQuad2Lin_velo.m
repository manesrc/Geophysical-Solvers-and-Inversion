function u_lin = convertQuad2Lin_velo(velo_quad,XP,TP,T)
% convert velocity field from quadratic to linear elements
% INPUT:
%   velo_quad: velocity field defined on quadratic elements
%   XP: coordinates of the linear mesh
%   TP: connectivity matrix of the linear mesh
%   T: connectivity matrix of the quadratic mesh
% OUTPUT:
%   u_lin: velocity field defined on linear elements


num_elem = size(T,1);
ndof_xlin = size(XP,1);
u_lin = zeros(ndof_xlin,2);

for j = 1:num_elem
    Te_quad = T(j,:);
    u_quad = velo_quad(Te_quad,:);
    Te_lin = TP(j,:);
    u_lin(Te_lin,:) = u_quad(1:length(Te_lin),:);
end 

end 