function [Ke,Ge,fe,mu22,temp22,rho22,Xee22] = makeElementMatrix(Xe, ...
   nVelocityDofPerElement,pospg,pespg,N,Nxi,Neta, ...
   nPressureNodesPerElement,NP, material,temp,Elem_LS,tolerance)
% input:
% Xe                        coordinates of the element nodes
% nVelocityDofPerElement    number of velocity dof per element
% pospg                     position of the gauss points
% pespg                     weights of the gauss points
% N, Nxi, Neta               shape functions and derivatives at gauss points
% nPressureNodesPerElement  number of pressure dof per element
% NP                       shape functions for pressure at gauss points
% material                  structure containing material properties
% temp                      nodal temperature values
% Elem_LS                   level set value for the element (if any)
% tolerance                 tolerance for level set evaluation
% output:
% Ke                        element stiffness matrix
% Ge                        element divergence matrix
% fe                        element force vector
% temp22, rho22, mu22      values of temperature, density, and viscosity at gauss points 
%                           (cross-validating computation of properties)
% Xee22                     positions where the above values have been computed (cross-validating mesh)


% initialisation
Ke = zeros(nVelocityDofPerElement,nVelocityDofPerElement);
Ge = zeros(nPressureNodesPerElement,nVelocityDofPerElement);
fe = zeros(nVelocityDofPerElement,1);

% number of integration points
ngaus = size(pospg,1);

% material properties
rho_ref = material.rho_ref; 
mu_ref  = material.mu_ref; 
g   = material.gravity;
alpha = material.alpha; 
beta = material.beta; 
T_0 = material.T_ref0;
mu_max = material.mu_max; 
A_d = material.disloc; 
eps_II = material.second_inv; 
flow_law = material.flow_law; 
energy  =material.energy; 
volume = material.volume;
gas_const = material.gas_constant; 
L_ref = material.L_ref;

mu22 = [];
rho22 = [];
temp22 = [];
Xee22 = [];


% loop in gauss points
for igaus = 1:ngaus   
   depth = (1-N(igaus,:) * Xe(:,2)) * L_ref; % [m]
   if size(temp,2) == 2
       LS_ig = N(igaus,:) * Elem_LS;
       if LS_ig < -tolerance
           T_1 = N(igaus,:) * temp(:,2);        % belongs to Omega_2
           mu_dim = material.mu2;
       elseif LS_ig > tolerance
           T_1 = N(igaus,:) * temp(:,1) ;       % belongs to Omega_1
           mu_dim = material.mu1;   
       else
           T_1 = 1523; % [K]
           mu_dim = material.mu2;
       end
   elseif Elem_LS ==  1
       mu_dim = material.mu1;
       T_1 = N(igaus,:)*temp;       
   elseif Elem_LS == 2
       mu_dim = material.mu2;
       T_1 = N(igaus,:)*temp;
       
   end
   press_lin = rho_ref*abs(material.gravity_units)*depth*(1e-6); % [MPa]  = [kg/m^3]*[m/s^2]*[m] * (1e-6 MPa / 1 Pa)
   rho_dim = rho_ref*(1- alpha*(T_1 - T_0) + beta*(press_lin) );        % [kg/m^3]
   rho_adim = rho_dim/rho_ref; 
   %calc_mu = A_d^(-1/flow_law)*eps_II^((1/flow_law)-1)*exp( (energy + press_lin*volume) / (flow_law*gas_const*T_1) ); % [MPa s]
   %mu_dim = min(mu_max,(1e3)*calc_mu); % [Pa s]
   mu_adim = mu_dim/mu_ref;

   mu22 = [mu22; mu_dim];
   rho22 = [rho22; rho_dim];
   temp22 = [temp22; T_1];
   Xee22 = [Xee22; N(igaus,:)*Xe];

   jacob = [Nxi(igaus,:) *Xe(:,1)  Nxi(igaus,:) *Xe(:,2); ...
            Neta(igaus,:)*Xe(:,1)  Neta(igaus,:)*Xe(:,2)];
   dvolu = pespg(igaus)*det(jacob);
   
   % shape functions derivatives in cartesian coordinates
   res = jacob\[Nxi(igaus,:);Neta(igaus,:)];   
   
   % resolution in terms of the stress-strain form ∫ eps * C * eps (page 278 
    % Donea-Huerta)
   B = [reshape([1;0]*res(1,:),1,nVelocityDofPerElement); ...
       reshape([0;1]*res(2,:),1,nVelocityDofPerElement);...
       reshape([res(2,:);res(1,:)],1,nVelocityDofPerElement)];
   C = eye(3);
   C(1:2,1:2) = 2*C(1:2,1:2);

   % diffusion matrix:
   %Ke = Ke + mu * (Nx'*Nx + Ny'*Ny) * dvolu;
   Ke = Ke + mu_adim*(B'*(C*B)) * dvolu;

   % divergence matrix:
   % divergence of shape functions
   dN = reshape(res,1,nVelocityDofPerElement);
   % element matrices
   Ge = Ge - NP(igaus,:)' * dN * dvolu;
   
   N_igaus = reshape([0;1]*N(igaus,:),1,nVelocityDofPerElement);
   fe = fe + N_igaus'*(rho_adim*g*dvolu);
end


%    % gradient of the shape functions (cartesian coords)
%    Nx = [reshape([1;0]*res(1,:),1,nVelocityDofPerElement); ...
%          reshape([0;1]*res(1,:),1,nVelocityDofPerElement)];
%    Ny = [reshape([1;0]*res(2,:),1,nVelocityDofPerElement); ...
%          reshape([0;1]*res(2,:),1,nVelocityDofPerElement)];