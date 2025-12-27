function matrix_M = compute_diff_meshes(DOF_fine,DOF_coarse,InfoMesh)
% Matrix M is a matrix to extend to DOF_G the originally q flux
% parametrized with ndofQ DOF
% INPUTS: 
% DOF_G: number of in the lower boundary of the mesh
% ndofQ: number of independent parameters describing q
% InfoMesh: Mesh information, 
% InfoProblem: Problem information.
%%%%%%%%%%%%%%%%%%%%%%%% CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_fine = length(DOF_fine);
num_coarse = length(DOF_coarse);
matrix_M = zeros(num_fine,num_coarse);
iniX = InfoMesh.ini_x; 
finX = InfoMesh.fin_x;
elemType = InfoMesh.elemType_LAB;   % line type of element (1-D)
% compute the integration points
num_gp = 2; 
num_nodes = 2;      % using linear elements 
Xfine = linspace(iniX,finX,num_fine);      % fine discretization [the size of the g vector: integration of q]
Xcoarse = linspace(iniX,finX,num_coarse);      % coarse discretization [the size of the q vector]

[chigp_fine, wgp_fine] = quadrature(elemType,num_gp);    % chigp_fine belongs to finer discr
[N_fine,~,~] = shapeFunctions(elemType,num_nodes,chigp_fine);        

for ii = 1:num_coarse-1           % coarse mesh nodes: columns
    indices_coarse = [ii ii+1];
    elem_coarse = Xcoarse(indices_coarse);     % find element in coarse mesh 
    % these are the nodes in the fine mesh between elem_coarse(1) and elem_coarse(2)
    relev_nodes_fine = DOF_fine((Xfine >= min(elem_coarse) & Xfine <= max(elem_coarse)));
    for jj = 1:length(relev_nodes_fine)-1
        indices_fine = [relev_nodes_fine(jj) relev_nodes_fine(jj+1)]; % node to the right (fine mesh)
        X_elem_fine = Xfine(indices_fine);   % element coordinates fine mesh
        xgp_fine = N_fine*X_elem_fine';     % global X coordinates of Gauss points
        chigp_coarse = interp1(elem_coarse, [-1 1], xgp_fine); % find isoparametric xgp_fine in coarse discretization
        [N_coarse, ~, ~] = shapeFunctions(elemType, num_nodes, chigp_coarse); % shape functions in coarse disc  
        long_elem_fine = X_elem_fine(2)-X_elem_fine(1);
        for mm = 1:num_nodes
            for qq = 1:num_nodes
                N_ii = N_coarse(:,mm);
                N_jq = N_fine(:,qq); % values of shape functions for node jj(qq)
                D_ij = compute_integral(N_ii,N_jq,wgp_fine,long_elem_fine);
                matrix_M(indices_fine(qq),indices_coarse(mm)) = matrix_M(indices_fine(qq),indices_coarse(mm)) + D_ij;
                %disp(['Componente M (',num2str(indices_fine(qq)),',',num2str(indices_coarse(mm)),')=',num2str(matrix_M(indices_fine(qq),indices_coarse(mm)))])
            end 
        end 
    end
end

end 
%%%%%%%%%%%%%%%%%%%%%%%% INTEGRAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D_ij = compute_integral(Ni,Nj,wgp,le)
    num_gp = length(wgp);
    D_ij = 0; % component of the matrix
    % compute integral numerically
    for kk = 1:num_gp
        D_ij = D_ij + wgp(kk) * Ni(kk) * Nj(kk) * (le/sum(wgp));
    end
end 
