function U_basis = VelocityBasis(str_name, example_plots, InfoProblem, InfoMaterial, InfoMesh, use_decay_model)
% Velocity basis: Computes the velocity basis for a rectangular mesh to be storing each component
% to be then used in a reduced order model
% OUTPUTS:
%   U_basis (matrix): Matrix where each column corresponds to the velocity field for a polygon
%   in the mesh.
% INPUTS:
%   str_name: Name of the .mat file to save the velocity basis.
%   example_plots: Structure containing plotting options.
%   InfoProblem: Structure containing problem information.
%   InfoMaterial: Structure containing material properties.
%   InfoMesh: Structure containing mesh information.
%   use_decay_model (0/1): 
%       - 1: Activates the code path with viscosity decay settings and LevelSet_closed2 on InfoMesh.X.
%       - 0: Uses the simpler LevelSet_closed + convertQuad2Lin path.

n_v_basis = InfoMesh.nel_x * InfoMesh.nel_y;
Xv = InfoMesh.X_v;
U_basis = zeros(size(Xv,1)*2,n_v_basis);
contador = 0; 
x_un1 = unique(Xv(:,1));
x_un = x_un1(1:2:end);   % for quadratic mesh 
y_un1 = unique(Xv(:,2));
y_un = y_un1(1:2:end);   % for quadratic mesh 
use_decay_model = logical(use_decay_model);

for jj = 1:length(y_un)-1
    ini_y = y_un(jj);
    fin_y = y_un(jj+1);
    for ii = 1:length(x_un)-1
        ini_x = x_un(ii);
        fin_x = x_un(ii+1);
        
        % use polygon 
        area_vertices = [ini_x ini_y; fin_x ini_y; fin_x fin_y; ini_x fin_y; ini_x ini_y];
        example_plots.LAB = area_vertices;
        
        % --- Logic Differentiator with IF statement ---
        if use_decay_model
            % Code path for the "most complex" function (VelocityBasis2)
            % Includes viscosity decay settings and LevelSet_closed2

            % center of polygon's coordinates
            InfoMaterial.ElemCenter = [(fin_x + ini_x)/2 (fin_y + ini_y)/2] / InfoProblem.L_ref;
            % length of polygon's side
            InfoMaterial.cellSide = min((fin_x - ini_x)/2, (fin_y - ini_y)/2) / InfoProblem.L_ref;
            % update Level Set
            InfoMesh.LS_new = LevelSet_closed(InfoMesh.X,area_vertices,example_plots);
            InfoMaterial.velo_basis = 1; % this flag enables the computation of decaying viscosity
        else
            % Code path for the simple VelocityBasis 
            % Uses LevelSet_closed and convertQuad2Lin
            % update Level Set
            LS_quad_new = LevelSet_closed(Xv,area_vertices,example_plots);
            InfoMesh.LS_new = convertQuad2Lin(LS_quad_new,InfoMesh.X,InfoMesh.T, InfoMesh.T_v);
            
            % NOTE: The simple path did NOT include InfoMaterial updates, so they are skipped here.
        end

        % LAB data
        InfoLAB.LABx = area_vertices(1:end,1);
        InfoLAB.LABy = area_vertices(1:end,2);
        
        % comput Stokes problem
        example_plots.fig_velo = 1;
        example_plots.LAB = [InfoLAB.LABx InfoLAB.LABy] / InfoProblem.L_ref;
        Temp_fake = [];
        p_fake = []; 
        
        [u_mant,~] = ComputeStokesProblem(Temp_fake, p_fake, InfoMesh,InfoProblem,InfoMaterial,example_plots);
        
        contador = contador + 1;
        % save in the basis
        U_basis(:,contador) = u_mant(:); 
    end 
end 

save(str_name,'U_basis'); %,'Writable',isWritable)
end