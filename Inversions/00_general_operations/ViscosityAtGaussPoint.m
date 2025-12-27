function mu = ViscosityAtGaussPoint(LS, Temp, pres, material)
% this function computes the viscosity at a Gauss point based on level set,
% temperature, pressure and material properties.
% INPUT:
%   LS        : level set value at the Gauss point
%   Temp      : temperature at the Gauss point (K)
%   pres      : pressure at the Gauss point (MPa)
%   material  : structure containing material properties and description
% OUTPUT:
%   mu        : viscosity at the Gauss point (Pa s)

% material description definition:
combinations_available = [material.mu_discont == 1 ...
    material.mu_pres_temp == 1 ...
    material.mu_smooth == 1];
material_definition = combinations_available * [1; 2; 3]; 
assert(isscalar(material_definition))

% material properties
mu1 = material.mu1; % [-]
mu2 = material.mu2; % [-]

% Level set;
dof_lithosphere = LS >= 0;
dof_asthenosphere = LS < 0;

if material_definition == 3
    error('Not implemented due to Temperature LS dependent')
elseif material_definition == 1
    mu = mu1 * dof_lithosphere + mu2 * dof_asthenosphere; 
elseif material_definition == 2
    mu = mu_const_pres_temp(pres,Temp,material);
else 
    error('The material description/combination is not implemented')
end 

end 

% density depending on temp and pressure
function mu_temp_pres_adim = mu_const_pres_temp(pressureMesh,TempMesh,InfoMaterial)
A_D = InfoMaterial.A_D; % MPa^{-n} s^{-1}
e_II = InfoMaterial.e_II;   % s^{-1}
E = InfoMaterial.E;     % J mol^{-1}
V = InfoMaterial.V;     % J MPa^{-1} mol^{-1}
n = InfoMaterial.n;     % [-]
R = InfoMaterial.R;     % J mol^{-1} K^{-1}
mu_max = InfoMaterial.mu_max;   % Pa s
exponential1 = (E+V*pressureMesh)./(n*R*TempMesh);
pre_multipl = (A_D^(-1/n)) * (e_II^((1/n)-1));
mu_temp_pres_adim = min(mu_max,pre_multipl*exp(exponential1)*1e6);  % 1e6 converts MPa s to Pa s

end 