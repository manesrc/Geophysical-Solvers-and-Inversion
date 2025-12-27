function [Ke,Ge,fe] = makeElementMatrix(Xe, nVelocityDofPerElement,pospg,pespg,N,Nxi,Neta, ...
   nPressureNodesPerElement,NP, g,rho_adim,mu_adim)
%   makeElementMatrix: builds the element matrices for the Stokes problem in 2D
%   OUTPUTS:
%   Ke : element diffusion matrix [nVelocityDofPerElement x nVelocityDofPerElement]
%   Ge : element divergence matrix [nPressureNodesPerElement x nVelocityD
%   fe : element force vector [nVelocityDofPerElement x 1]
%   INPUTS:
%   Xe : coordinates of the element nodes [nNodesPerElement x 2]
%   nVelocityDofPerElement : number of velocity DOF per element
%   pospg : positions of the integration points [nGausspoints x 2
%   pespg : weights of the integration points [nGausspoints x 1]
%   N : shape functions at the integration points [nGausspoints x nNodes
%   Nxi : derivatives of shape functions respect xi at the integration points
%         [nGausspoints x nNodes]
%   Neta : derivatives of shape functions respect eta at the integration points
%          [nGausspoints x nNodes]
%   nPressureNodesPerElement : number of pressure DOF per element
%   NP : shape functions for pressure at the integration points
%   g : gravity (adim)
%   rho_adim : density at the integration points (adim) [nGausspoints x 1]
%   mu_adim : viscosity at the integration points (adim) [nGausspoints x 1]

% initialisation
Ke = zeros(nVelocityDofPerElement,nVelocityDofPerElement);
Ge = zeros(nPressureNodesPerElement,nVelocityDofPerElement);
fe = zeros(nVelocityDofPerElement,1);

% number of integration points
ngaus = size(pospg,1);

% material properties
%g   = material.gravity; % adim


% loop in gauss points
for igaus = 1:ngaus   
    
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
   Ke = Ke + mu_adim(igaus)*(B'*(C*B)) * dvolu;

   % divergence matrix: divergence of shape functions
   dN = reshape(res,1,nVelocityDofPerElement);
   % element matrices
   Ge = Ge - NP(igaus,:)' * dN * dvolu;
   
   N_igaus = reshape([0;1]*N(igaus,:),1,nVelocityDofPerElement);
   fe = fe + N_igaus'*(rho_adim(igaus)*g*dvolu);
end