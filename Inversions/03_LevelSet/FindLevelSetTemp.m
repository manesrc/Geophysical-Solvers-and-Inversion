function LS_temp_Xv = FindLevelSetTemp2(LSTemp,nnodes,T,TP)
% OUTPUT: LS_temp_Xv is the value of the level set for the nodes in the quadratic mesh
% INPUTS:
% LSTemp : existing LS in the linear mesh
% nnodes: number of nodes in the quad mesh 
% TP : element connectivities in the linear mesh
% T : element connectivities in the quadratic mesh

LS_temp_Xv = zeros(nnodes,1);

assert(size(T,1) == size(TP,1))
assert(T(1,1) == TP(1,1))
num_elem = size(T,1);

zgp = []

for j = 1:num_elem
    Te_lin = TP(j,:);
    
    u_quad = quad_field(Te_quad,:);
    Te_lin = TP(j,:);
    lin_field(Te_lin,:) = u_quad(1:length(Te_lin),:);
end 

end 