function InfoMesh = generate_mesh_related_data(nel_x, nel_y, Long_X,maxDepth,InfoProblem)
% generate_mesh_related_data - Generates mesh related data for the simulation.
% Inputs:
%    nel_x     - Number of elements in the x-direction
%    nel_y     - Number of elements in the y-direction
%    Long_X    - Length of the domain in the x-direction [m]
%    maxDepth  - Maximum depth of the domain [m]
%    InfoProblem - Structure containing problem information, including L_ref to be in [m] as well
% Outputs:
%    InfoMesh  - Structure containing mesh related information

    % Generate Fixed mesh
    InfoMesh.ini_x = 0; 
    InfoMesh.fin_x = Long_X / InfoProblem.L_ref; 
    InfoMesh.ini_y = 0; 
    InfoMesh.fin_y = maxDepth/InfoProblem.L_ref; 
    % elements
    InfoMesh.nel_x = nel_x;
    InfoMesh.nel_y = nel_y; 
    InfoMesh.elemType = 1;
    % temperature mesh
    InfoMesh.nne = 4;
    addpath '02_mesh_cond_SF'
    [InfoMesh.X,InfoMesh.T] = CreaMalla_rectangulo(InfoMesh.nne, InfoProblem.L_ref, InfoMesh);
    % velocity mesh
    InfoMesh.ngp_pois=9;
    InfoMesh.nnv = 9;
    [InfoMesh.X_v,InfoMesh.T_v] = CreaMalla_rectangulo(InfoMesh.nnv, InfoProblem.L_ref, InfoMesh);
end