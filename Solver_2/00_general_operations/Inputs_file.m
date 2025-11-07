function [InfoProblem, InfoMaterial] = Inputs_file(InfoProblem)
    % Inputs_file: Defines problem and material properties for the simulation 
    % created to ease the problem definition at main
    % outputs:
    %   InfoProblem  : Structure containing problem properties
    %   InfoMaterial : Structure containing material properties
    % inputs:
    %   InfoProblem  : Structure containing problem properties to be updated

% Conditions and material properties
InfoProblem.T_sup = 293;    % [K]
InfoProblem.T_LAB = 1573;   % [K]
InfoProblem.k1 = 1; % [W/(m K)]
InfoProblem.s1 = 0; % [W/m^3] 
InfoProblem.grad_apprT2 = 0.5/1000; % [K/m]
InfoProblem.s2  = 0;    % s2 = rho2*f2 [W/m^3]
% for velocity purposes:
InfoMaterial.gravity_units = -9.81; % [m/s^2]
InfoMaterial.calorific = 1200;   % [J/(kg*K)]
InfoMaterial.mu1 = 1e24; % [Pa*s] 
InfoMaterial.mu2 = 1e20; % [Pa*s]
InfoMaterial.rho_ref = 3300; 
InfoMaterial.mu_ref = max(InfoMaterial.mu1,InfoMaterial.mu2);
% rho and mu definition: (scatter, smooth, temp-pres dependent)
InfoMaterial.rho_pres_temp = 1; InfoMaterial.rho_discont = 0; %,1
InfoMaterial.mu_discont = 1;  InfoMaterial.mu_pres_temp = 0; %,1
InfoMaterial.mu_smooth = 0;     % if mu_smooth == 1 define InfoMaterial.thickness (thickness in which mu_2 varies into mu_1)
% material definition subjected to mu_smooth == 1
InfoMaterial.thickness = 50*1000; % [m]
% material definitions for alpha and beta constants related to density function
InfoMaterial.alpha_dens = 1e-5; % [1/K]
InfoMaterial.beta_dens = 1e-5; % [1/MPa]
% Define k_variable with temp and pressure
InfoMaterial.k_var = 0;
% Info Material for k(depth) alpha (T) [Fei 95]
InfoMaterial.alphacoeff1 = 2.65032e-5;
InfoMaterial.alphacoeff2 = 9.19917e-9;
InfoMaterial.alphacoeff3 = -0.2712774;
% with this rho depends on alpha that is a function of the temperature 
% therefore alpha is not constant
InfoMaterial.rho_alphaT = 0; 

% k_rad (T) [Hoeffmeister 99]
InfoMaterial.k_rad1 = 1.753e-2;
InfoMaterial.k_rad2 = -1.0365e-4;
InfoMaterial.k_rad3 = 2.2451e-7;
InfoMaterial.k_rad4 = -3.407e-11;
% k_gen (depth) [Afonso 08]
InfoMaterial.k_0 = 2.5;     % [W/(m*K)]
InfoMaterial.a = 0.45;
InfoMaterial.gamma = 1.25;
InfoMaterial.K0 = 120; % [GPa]
InfoMaterial.derK0 = 4.5; % [GPa]

% Data for temperature and pressure dependent viscosity
InfoMaterial.A_D =1.1e5;  % [MPa^(-n)/s]
InfoMaterial.e_II = 1e-15;   % [1/s]
InfoMaterial.E = 5.3e5; % [J/mol]
InfoMaterial.V = 14; % [J / (MPa*mol)]
InfoMaterial.n = 3.5;    % [-]
InfoMaterial.R = 8.314; % [J/(mol*K)]
InfoMaterial.mu_max = 1e24;  % [Pa s]


% material description definition:
combinations_available = [(InfoMaterial.rho_pres_temp == 1 && InfoMaterial.mu_discont == 1) ...
    (InfoMaterial.rho_pres_temp == 1 && InfoMaterial.mu_pres_temp == 1) ...
    (InfoMaterial.rho_pres_temp == 1 && InfoMaterial.mu_smooth == 1) ...
    (InfoMaterial.rho_discont == 1 && InfoMaterial.mu_discont == 1)];
Mat_defin = combinations_available * [1; 2; 3; 4]; 

str1 = 'rho=f(T,p) & mu discontinuous';
str2 = 'rho=f(T,p) & mu=g(T,p)';
str3 = 'rho=f(T,p) & mu smooth transition from mu_1 to mu_2 (use a transition thickness parameter)';
str4 = 'rho and mu discontinuous';

if Mat_defin == 1
    str = str1;
elseif Mat_defin == 2
    str = str2;
elseif Mat_defin == 3
    str = str3;
else
    str = str4; 
end 

disp('Summary of properties: ')
disp(['Tsup = ', num2str(InfoProblem.T_sup),', TLAB = ', num2str(InfoProblem.T_LAB),', s1 = s2 = ', num2str(InfoProblem.s1),', mu_1 = ',num2str(InfoMaterial.mu1),', mu_2 = ',num2str(InfoMaterial.mu2)])
if InfoMaterial.k_var == 1
    str10 = ['conductivity k_1 = f(p,T)'];
else
    str10 = ['conductivity k_1 = ', num2str(InfoProblem.k1)];
end 
disp(str10)
disp(['material description (',num2str(Mat_defin),'): ', str])

end 