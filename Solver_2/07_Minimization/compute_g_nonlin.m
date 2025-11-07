function g_nl = compute_g_nonlin(thickness, InfoMesh)
% OUTPUT: velocity modified
% INPUT: delta (linear combination for Ubasis)
%       thickness (thickness in which velocity is attenuated)
%       Ubasis velocity basis
%       InfoMesh many mesh-related data.. mesh connectivities, etc
% explanation: this function computes the velocity with Ubasis and alpha and then erases 
% the DOFs outside Omega2 and applies a gradual attenuation close to LAB
% lists are formed with mesh "X" ; velocities goes with mesh "X_v"
    X_v = InfoMesh.X_v;
    g_nl = zeros(size(X_v)); 
    ndof_u_mesh = size(X_v,1);
    attenuation_function = zeros(ndof_u_mesh,1);
    T_v = InfoMesh.T_v;    
    % temp mesh
    T = InfoMesh.T; % X = InfoMesh.X;
    % num of elements
    num_elements = size(T,1);
    % Level-Set re-scaled to distance    
    LS_distance = InfoMesh.LS_val; 
    % adjust a curve to modify velocities close to LAB
    factor1 = 1; 
    chi = 5 / (thickness * factor1); 
    isop_points_velo_mesh = [-1 -1; 1 -1; 1 1; -1 1; 0 -1; 1 0; 0 1; -1 0; 0 0];
    N_velo_mesh = shapeFunctions(1,size(T,2),isop_points_velo_mesh);
    for ii = 1:num_elements
        Te = T(ii,:); % connectivity in Temp mesh
        LS_distance_e = LS_distance(Te); % Level-Set values
        LS_velo_mesh = N_velo_mesh * LS_distance_e; % convert LS to velo mesh
        Te_v = T_v(ii,:); % connectivity in velo mesh
        if sum(sign(LS_velo_mesh)) == length(LS_velo_mesh) % all nodes in Omega1
            %attenuation_function(Te_v) = 1e-2;
            factor123 = 1/400;
            alpha2 = 5/thickness; 
            %attenuation_function(Te_v) = factor123 - factor123 * exp(-alpha2 * LS_velo_mesh);
            attenuation_function(Te_v) = 0; 
        else
            Te_inf = Te_v(LS_velo_mesh < 0); 
            bool_LS = LS_velo_mesh < 0;
            attenuation_function(Te_inf) = 1 - exp(chi*LS_velo_mesh(bool_LS));
        end 
    end 
    % velocity of interest  
    g_nl(:,1) = attenuation_function;
    g_nl(:,2) = attenuation_function;
    
end 