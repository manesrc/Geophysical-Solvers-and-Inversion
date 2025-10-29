function [Ge1,Me2,be1,me2,Me2_int] = EleMatNitsche(Xe,InfoProblem,InfoMesh,k,problem_int,Ar_phys,intersection,normal2elem,plot_points_interface)
% Outputs:
% Ge1: Elemental gradient matrix Ge1 = [ ∫ (∂N/∂x * n) * N d Gamma]
% Me2: Elemental mass matrix Me2 =  [ß ∫ N * N d Gamma]
% be1: Elemental gradient vector be1 = [ ∫ (∂N/∂x * n) * T_LAB d Gamma ]
% me2: Elemental mass vector me2 = [ ∫ N * T_LAB d Gamma ]
% Me2_int: Elemental mass matrix without considering ß, Me2_int = [ ∫ N * N d Gamma]
% Inputs: 
% Xe: coordinates of the element nodes
% InfoProblem: data of the problem as Length, Temperatures, material parameters, etc
% InfoLAB: Information on the LAB position
% InfoMesh: X matrix, T matrix, number of elements, etc
% k: material parameter of the problem of interest
% problem_int: Omega 1 ==  1, Omega 2 == 2
% Ar_phys: part of the element in problem_int
% intersection: point of intersection element-interface, P1 and P2 in isoparametric coordinates [1 x 4]
% normal2elem: Normal to the element pointing outwards Omega1

% Allocation
ndof_e = size(Xe,1);           
Ge1 = zeros(ndof_e);
Me2_int = zeros(ndof_e);
be1 = zeros(ndof_e,1);
me2_int = zeros(ndof_e,1);

% Number of GP integrate in LAB interface line
NumberOfPointsLAB = InfoMesh.ngp_interface; 
% intersection coordinates
in_coord = intersection(1:2);
out_coord = intersection(3:4);
% compute the length of the intersection line
isop_bound = [-1 1];
% X,Y boundaries
x_bound = [min(Xe(:,1)) max(Xe(:,1))];
y_bound = [min(Xe(:,2)) max(Xe(:,2))];
% intersection points in X,Y coord. to calculate length of intersection
chi2x = interp1(isop_bound,x_bound,[in_coord(1) out_coord(1)]);
eta2y = interp1(isop_bound,y_bound,[in_coord(2) out_coord(2)]);
% length_intersection to then calculate ß stabilization parameter
length_intersection = sqrt( (chi2x(1) - chi2x(2))^2 + (eta2y(1) - eta2y(2))^2 );
% compute the position of the Gauss Points in a line [interface] on 2D space
[chieta_gp,wgp] = quad2D_line(NumberOfPointsLAB,in_coord,out_coord);

% Boundary condition to impose weakly
Temp_LAB = InfoProblem.T_LAB;

%% compute the shape functions to integrate
elemType2D = InfoMesh.elemType; % 2D elem type
NofNodes2D = InfoMesh.nne;          % 2D elem type
% calculate shape functions in chieta_gp
[N,Nxi,Neta] = shapeFunctions(elemType2D,NofNodes2D,chieta_gp);

% plot the chieta_gp in figure 120 [graphical check]
if plot_points_interface == 1
    figure(120); hold on; scatter(N*Xe(:,1),N*Xe(:,2),3,'black','x'); %set(h,'AutoUpdate','off')
end 

%% Loop in Gauss points
% to compute  Nitsche factor[beta]
ngaus = size(wgp,1);

for igaus = 1:ngaus
    % isop. shape func. and derivatives
    N_igaus = N(igaus,:);
    Nxi_igaus = Nxi(igaus,:);
    Neta_igaus = Neta(igaus,:);
    % Jacobian 
    Jacob = [Nxi_igaus;Neta_igaus]*Xe;
    dN_iso = [Nxi_igaus;Neta_igaus];
    % gradient in terms of X,Y
    dN_dxg = Jacob\dN_iso;
    % surface of integration
    dline = (length_intersection/sum(wgp))*wgp(igaus);          
    %dline = (length_intersection/2)*wgp(igaus);      %  esto estaría mal     
    % normal to the surface of integration
    if problem_int == 1
        normal = normal2elem';
    else
        normal = -normal2elem';
    end 
    % Ge1 = int_\Gamma [ k * N_i * ( Nabla N_j * normal ) +  N_j * (k* Nabla N_i * normal ) ] d (\Gamma)
    Ge1 = Ge1 + k * ((N_igaus') * (dN_dxg'*normal)' + (dN_dxg'*normal) * N_igaus) * dline;  
    % Ke2 = \int_\Gamma N_i * (N_j)
    Me2_int = Me2_int + (N_igaus'*N_igaus)*dline;
    % fe_1 = int_\Gamma u_d * ( Nabla N_j * normal ) d (\Gamma)
    be1 = be1 + k * Temp_LAB * (dN_dxg'*normal) * dline;
    % fe_2 = int_\Gamma (u_d * N_j) d (\Gamma)   
    me2_int = me2_int + Temp_LAB * N_igaus' * dline;      
end

% Dolbow & Harari: 
el_size = (1/InfoMesh.nel_x)*(1/InfoMesh.nel_y); % units [0,1] x [0,1]
ele_domain_int = el_size*Ar_phys;           % units [0,1] x [0,1], Ar_phys is unitless

if ele_domain_int == 0
    beta = 0;
    keyboard
else 
    d = size(Xe,2);         % dimension of the problem (2D)
    if (size(N_igaus,2) == 9) || (size(N_igaus,2)==7)
        p = 2; % interpolation degree
    else
        p = 1; 
    end 
    higher_degree_interp = InfoProblem.c*(2* p *(p-1+d))/d;
    beta =  higher_degree_interp * k * length_intersection / ele_domain_int; 
end

%beta = alpha * (k/(Ar_phys*h));
Me2 = beta*Me2_int;
me2 = beta*me2_int;

end
