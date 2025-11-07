function mu = ViscosityAtGaussPoint2(LS, Temp, pres, material)
% ViscosityAtGaussPoint: Computes the viscosity at Gauss points based on material 
% properties and distance to Level-Set (LS = 0 is the LAB)
% outputs:
%   mu      : Viscosity at Gauss points
% inputs:
%   LS       : Level-Set function values
%   Temp     : Temperature at Gauss points
%   pres     : Pressure at Gauss points
%   material : Material properties structure

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
LS_at_lithosphere = LS >= 0;
LS_at_asthenosphere = LS < 0;

if material_definition == 3
    error('Not implemented due to Temperature LS dependent')
elseif material_definition == 1
    mu = mu1 * LS_at_lithosphere + mu2 * LS_at_asthenosphere; 
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