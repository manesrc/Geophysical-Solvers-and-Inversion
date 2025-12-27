function lin_field = convertQuad2Lin(quad_field,XP,TP,T)
% convertQuad2Lin: converts results from a quadratic mesh to a linear mesh
% OUTPUT: 
% u_lin velocity results in the linear mesh
% INPUTS:
% quad_field : existing results in the quadratic mesh at points X (not used in the code)
% XP : nodes in the linear mesh
% TP : element connectivities in the linear mesh
% T : element connectivities in the quadratic mesh

assert(size(T,1) == size(TP,1))
assert(T(1,1) == TP(1,1))

num_elem = size(T,1);
ndof_xlin = size(XP,1);
lin_field = zeros(ndof_xlin,size(quad_field,2));

for j = 1:num_elem
    Te_quad = T(j,:);
    u_quad = quad_field(Te_quad,:);
    Te_lin = TP(j,:);
    lin_field(Te_lin,:) = u_quad(1:length(Te_lin),:);
end 

end 