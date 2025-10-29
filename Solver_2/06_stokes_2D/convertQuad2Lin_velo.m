function u_lin = convertQuad2Lin_velo(velo_quad,XP,TP,T)

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