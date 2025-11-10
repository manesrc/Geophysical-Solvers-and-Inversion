function [g_vect,line_int_vect] = compute_g_mean(DOF_G,InfoMesh,InfoProblem)
% compute the mean value of g over the interface LAB
% input:
% DOF_G: vector with the DOFs associated to the interface LAB
% InfoMesh: structure with the mesh information
% InfoProblem: structure with the problem information
% output:
% g_vect: vector with the mean value of fluxes over the bottom boundary
% line_int_vect: vector with the integral of g over the interface LAB
numDOF_inf = size(DOF_G,1);
x_length_elem = (1/(numDOF_inf-1)); % element size
% 1D shape functions
T = InfoMesh.T;
nne_line = sqrt(size(T,2));
[xgp,wgp] = quadrature(InfoMesh.elemType_LAB,nne_line);
[N,~,~] = shapeFunctions(InfoMesh.elemType_LAB,nne_line,xgp);

line_int_vect = zeros(numDOF_inf,1);
for jj = 1:(numDOF_inf-1)
    %%% si se corren casos con elementos cuadráticos no va a funcionar
    if nne_line == 3
        Te = [2*jj-1 2*jj 2*jj+1];
    else
        Te = [jj jj+1];
    end 
    g_e = zeros(size(N,2),1);
    for kk = 1:size(wgp,2)
        g_e = g_e + N(kk,:)'*wgp(kk) * (x_length_elem/2);
    end
    line_int_vect(Te,1) = line_int_vect(Te,1) + g_e;
end 
g_vect = InfoProblem.q2*line_int_vect;
end 