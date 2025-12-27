function [Temp, pressure] = estimate_Temp_pres(X, nel_x, nel_y, LS_geom, InfoMaterial, InfoProblem)
    % estimate_Temp_pres: Estimate temperature and pressure field based on
    % level set geometry and material/problem information.
    % INPUT:
    %   X             : Nodal coordinates
    %   nel_x         : Number of elements in x-direction
    %   nel_y         : Number of elements in y-direction
    %   LS_geom       : Level set values at nodal points
    %   InfoMaterial  : Structure containing material properties
    %   InfoProblem   : Structure containing problem settings   

    % OUTPUT:
    %   Temp          : Estimated temperature field at nodal points
    %   pressure      : Estimated pressure field at nodal points
    
    % initialize
    Temp = zeros(size(X,1),1);
    rho0 = InfoMaterial.rho_ref;
    g = abs(InfoMaterial.gravity_units);
    % calculate pressure
    pressure = rho0 * g * (max(X(:,2)) - X(:,2) ) * (1e-6); % MPa
    % use available temperatures
    T_LAB = InfoProblem.T_LAB;
    T_sup = InfoProblem.T_sup;
    T_inf = InfoProblem.T_inf;

    tot_dof = size(X,1);

    ind_nodes = 1:1:tot_dof;
    for i = 1:(nel_x+1)
        dof_i = ind_nodes(i);
        LS_inf = LS_geom(dof_i);
        nodes_col_i = dof_i:(nel_x+1):tot_dof;
        LS_sup = LS_geom(nodes_col_i(end));
        for j = 1:length(nodes_col_i)
            LS_ij = LS_geom(nodes_col_i(j));
            if LS_ij < 0
                slope = (T_inf-T_LAB) / LS_inf;
            else
                slope = (T_sup - T_LAB) / LS_sup;
            end 
            Temp( nodes_col_i(j) ) = T_LAB + slope * LS_geom( nodes_col_i(j) );
        end 
    end

end         