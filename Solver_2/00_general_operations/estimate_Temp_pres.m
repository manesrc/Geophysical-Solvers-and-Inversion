function [Temp, pressure] = estimate_Temp_pres(X, nel_x, nel_y, LS_geom, InfoMaterial, InfoProblem)
% estimate_Temp_pres: Estimates temperature and pressure at nodes based on Level-Set 
% geometry and material/problem info
% outputs:
%   Temp     : Estimated temperature at nodes (bi-linear gradient from LS, note LS = 0 is the LAB)
%   pressure : Estimated pressure at nodes
% inputs:
%   X            : Coordinates of nodes
%   nel_x        : Number of elements in x-direction
%   nel_y        : Number of elements in y-direction
%   LS_geom      : Level-Set function values at nodes
%   InfoMaterial : Material properties structure
%   InfoProblem  : Problem properties structure

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

% function [Temp,p] = estimate_Temp_pres(X,nel_x,nel_y,LS_geom,InfoMaterial,InfoProblem)
% % initialize
% Temp = zeros(size(X,1),1);
% rho0 = InfoMaterial.rho_ref;
% g = abs(InfoMaterial.gravity_units);
% % calculate pressure
% p = rho0 * g * (max(X(:,2)) - X(:,2) ) * (1e-6);% MPa
% % use available temperatures
% T_LAB = InfoProblem.T_LAB;
% T_sup = InfoProblem.T_sup;
% T_inf = InfoProblem.T_inf;
% 
% ind = 1:1:nel_y+1;
% for ii = 1:(nel_x+1)
%     nodes_col_ii = ii:nel_x+1:size(X,1);
%     LS_ii = LS_geom(nodes_col_ii);
%     if sum(LS_ii == 0) == 0
%         last_node_Omega2 = max(ind(LS_ii < 0));
%         first_node_Omega1 = min(ind(LS_ii > 0));
%         assert(last_node_Omega2+1 == first_node_Omega1)
%         % interpolate to find Gamma_LAB (at X=X(ii,2))
%         LS_var = [LS_ii(last_node_Omega2) LS_ii(first_node_Omega1)];
%         Y_var = [X(nodes_col_ii(last_node_Omega2),2) X(nodes_col_ii(first_node_Omega1),2)];
%         Y_zero = interp1(LS_var,Y_var,0);   % distance to bottom
%         node2mod = [last_node_Omega2 first_node_Omega1] * (abs(Y_var - Y_zero) == min(abs(Y_var - Y_zero)) )'; % modify the temperature of the closest node to LAB to ensure T=TLAB at yLAB
%     else
%         ind0 = ind(LS_ii == 0);
%         Y_zero = X(nodes_col_ii(ind0),2);
%         node2mod = [];
%     end 
% 
%     grad_inf = 0.5/1000;
%     grad_sup = (T_LAB - T_sup) / (max(X(:,2)) - Y_zero);
%     depth_LAB = max(X(:,2)) - Y_zero; 
%     for jj = 1:length(nodes_col_ii)
%         node_jj = nodes_col_ii(jj);
%         depth_jj = max(X(:,2)) - X(node_jj,2);
%         if depth_jj > depth_LAB % below LAB
%             Temp(node_jj) = T_LAB + grad_inf * (depth_jj - depth_LAB); 
%         else
%             Temp(node_jj) = T_sup + grad_sup * depth_jj;
%         end 
%         if jj == node2mod   % modify temperature on node closest to y_LAB
%             pend_int = grad_inf; 
%             assert(pend_int > 0)
%             Temp(node_jj) = T_LAB + pend_int * (depth_jj - depth_LAB); 
%         end 
%     end 
% end
% 
% end 