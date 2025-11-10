function [chi_pg,wpg,Ar_phys_el] = ModifyQuadrature_LevelSet(InfoMesh,problem_int,LSe,tolerance)
%% compute more gauss points to see which of them is inside the element
% and which outside depending on the level set function
% Inputs:
%   InfoMesh: structure containing mesh information
%   problem_int: 1 - Lithosphere, 2 - Asthenosphere
%   LSe: level set values at the element nodes
%   tolerance: tolerance to consider a point on the interface
% Outputs:
%   chi_pg: locations of the gauss points
%   wpg: modified weights of the gauss points
%   Ar_phys_el: percentage of the element that belongs to the physical domain

elemType = InfoMesh.elemType;                                           % quad
NumberOfGaussPoints = InfoMesh.ngp_Nits;                         % more gauss points than normal
[chi_pg,pespg] = quadrature(elemType,NumberOfGaussPoints); % obtain location of GP [chi,eta]
%NumberOfElementNodes = InfoMesh.nne;                            % number of element nodes
NumberOfElementNodes = 4;                            % Linearize shape functions
[N,~,~] = shapeFunctions(elemType,NumberOfElementNodes,chi_pg); % shape functions to calculate LS_{GP}

% set to zero those LS values smaller than the tolerance in abs value
for jj = 1:length(LSe)
    if abs(LSe(jj)) < tolerance
        LSe(jj) = 0;
    end 
end

LS_gp = N*LSe(1:NumberOfElementNodes);  % LSe are the LS values of the nodes and depend on (X,Y) 
wpg = pespg;
for ii=1:length(LS_gp)
    if (problem_int == 1) && LS_gp(ii) < -tolerance
        wpg(ii) = 0;
    elseif problem_int == 2 && LS_gp(ii) > tolerance
        wpg(ii) = 0;
    elseif abs(LS_gp(ii)) < tolerance       % the points belong in the interface (will appear twice)
        wpg(ii) = pespg(ii)/2;
    end 
end 
% compute the percentage of the element in the phys. domain
Ar_phys_el = sum(wpg)/sum(pespg); 

end