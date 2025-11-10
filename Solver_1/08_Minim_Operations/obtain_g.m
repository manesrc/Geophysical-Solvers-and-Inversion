function [g_perturbation,g_pert_sq] = obtain_g(B_0,v_0,method1,ratio_eig)
% obtain_g: function to obtain the perturbation g (flux at the bottom boundary) 
% that minimizes the flux across the interphase Gamma
%
% INPUT:
% B_0                 : matrix B_0 of the rectangular system B_0*q
% v_0                 : vector v_0 of the rectangular system B_0*q
% method1             : method to solve the system
%                       1: rectangular system
%                       2: square system
%                       4: both methods
% ratio_eig          : ratio between the biggest and smallest eigenvalue to
%                       consider in the minimization
% OUTPUT:
% g_perturbation     : perturbation g obtained with the rectangular system
% g_pert_sq          : perturbation g obtained with the square system (if
%                       method1 = 4)

%% main code
if method1 == 1 % solve the rectangular system 
    g_perturbation = compute_perturbation(B_0,v_0,ratio_eig);
    g_pert_sq = [];
elseif method1 == 2    % solve the square system
    [~,S12,~] = svd(B_0'*B_0);
    D_mat2 = eye(size(B_0,2));
    alpha = sqrt(S12(1,1))*ratio_eig;
    mat1_q = (B_0'*B_0 + (alpha^2)*D_mat2);
    vect1_q = (B_0'*v_0);
    g_perturbation = mat1_q\vect1_q ; % values q minimizing the flux across Gamma
    g_pert_sq = [];     % not really 
elseif method1 == 4     % run both method 1 and 2 and save results
    % first compute solution of the square system
    [~,S12,~] = svd(B_0'*B_0);
    D_mat2 = eye(size(B_0,2));
    alpha = sqrt(S12(1,1))*ratio_eig; 
    mat1_q = (B_0'*B_0 + (alpha^2)*D_mat2);
    vect1_q = B_0'*v_0;
    % square solution
    g_pert_sq = mat1_q\vect1_q;
    % rectangular solution 
    g_perturbation = compute_perturbation(B_0,v_0,ratio_eig);
end 

% solution 
normV0 = norm(v_0);
 disp(['Norm of [ v_0 ] = ',num2str(normV0)])
normRes = norm(v_0 - B_0*g_perturbation);
disp(['Norm of [ B_0*q - v_0 ] = ',num2str(normRes)])
disp(['Performance of the minimization 100* [ 1- (Norm of [v_0-B_0*g_solution]/Norm of [v_0]) ] = ',num2str(100*round(1-(normRes/normV0),3)),'%'])

end     

%% rectangular solution 
function pert = compute_perturbation(B_0,v_0,diff_betw_eigenval)
%v_0_fin = v_0 - B_0 * g_mean;   % compute the vector minus the initial value   
[U,S,V] = svd(B_0); % singular value decomposition of matrix

if size(S,2) == 1
    vect1 = 1;
else
    vect1 = diag(S)/S(1,1);
end

vect2 = vect1 > diff_betw_eigenval;
n_red = sum(vect2);

z1 = U'*v_0;        

y1 = zeros(n_red,1);

for ii = 1:n_red
    y1(ii) = z1(ii)/S(ii,ii);
end 
pert = V(:,1:n_red)*y1;
end 

