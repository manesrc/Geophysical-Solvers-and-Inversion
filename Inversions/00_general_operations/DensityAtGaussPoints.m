function rho = DensityAtGaussPoints(LS, Temp, pres,T_sup,material)
% this function computes the density at gauss points depending on the 
% material definition

% INPUT:
% LS: level-set values at gauss points
% Temp: temperature at gauss points
% pres: pressure at gauss points
% T_sup: surface temperature
% material: structure containing material properties

% OUTPUT:
% rho: density at gauss points

% material definition:
combinations_available = [material.rho_pres_temp == 1 ...
    material.rho_discont == 1];
material_definition = combinations_available * [1; 2]; 
assert(isscalar(material_definition))
    
% define rho: 
rho0 = material.rho_ref;

% Level-Set
dof_lithosphere = LS >= 0;
dof_asthenosphere = LS < 0;

if material_definition == 1 % rho: func(pres,T)
    % alpha coefficient:
    alpha_T = material.rho_alphaT; 
    beta = material.beta_dens;
    if alpha_T == 1
        a1 = material.alphacoeff1;
        b1 = material.alphacoeff2;
        c1 = material.alphacoeff3; 
        alpha1 = a1 + b1 * Temp + c1 * Temp.^-2;
        rho = rho0* (1 - alpha1 .* (Temp - T_sup) + beta * pres);
    else
        alpha = material.alpha_dens;
        rho = rho0* (1 - alpha * (Temp-T_sup) + beta * pres);
    end     
elseif material_definition == 2     % rho: {rho1, rho2}
    % position of gauss points
    rho1 = material.rho1; % [-]
    rho2 = material.rho2; % [-]
    rho = rho1*(dof_lithosphere)+rho2*(dof_asthenosphere);
else 
    error('The material description/combination is not implemented')
end

end 