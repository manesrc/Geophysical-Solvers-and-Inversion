function [U_basis2,alpha_min,velo_min,norm_dif,contador_svd,final_size] = SVD_ordered_Ubasis(u_ref,reduction,nel_x,nel_y,U_basis,InfoMesh,do_plots)
% OUTPUT: 
%       U_basis_new = modified velocity basis
%       alpha_min = linear combination of U_basis_new approximating u_ref
%       velo_min = U_basis_new * alpha_min
%       norm_dif = || u_ref - velo_min || / ||u_ref||
%       contador_svd = number of components in the basis disregarded
%       final_size = number of components in U_basis_new
% INPUT:
%       u_ref = velocity of reference or to approximate
%       reduction case (see below)
%       nel_x, nel_y: number of elmenets in the mesh 
%       U_basis: Original velocity basis to modify 
%       InfoMesh: characteristics of the mesh
%       doplot: 
% reduction cases: 
% 0 = no reduction original U_basis, 
% 0.5 trim last component of the basis, 
% 1 trim last column, 
% 2 trim first and last, 
% 3 use the symmetric base, 
% 4 use basis sumed every four elements in squares, 
% 5 use the linear independent SVD decomposition,
% 6 symmetric base minus "first and last columns" (technically only takes one column out and the other dissappears by simmetry); 
% 10 is taking out those elements above LAB; 
% 11 takes out those above LAB plus last column
% 16 use elements summed in 4 and take out last column
% 17 from base grouped by four take out those from Lithosphere

    if reduction == 0   % 
        U_basis2 = U_basis; 
    elseif reduction == 0.5 % trim the last component of the base
        U_basis2 = U_basis(:,1:end-1);
    elseif reduction == 1    % trim last column
        % trim last column 
        indices2trim1 = nel_x:nel_x:nel_x*nel_y;
        U_basis2 = U_basis;
        U_basis2(:,indices2trim1) = [];
    elseif reduction == 2    % trim first and last columns
        % trim last column 
        indices2trim1 = nel_x:nel_x:nel_x*nel_y;
        % trim first column
        indices2trim2 = 1:nel_x:nel_x*nel_y;
        indices2trim = [indices2trim1 indices2trim2];
        U_basis2 = U_basis;
        U_basis2(:,indices2trim) = [];
    elseif reduction ==  3
        U_basis2 = generate_symmetric_base(U_basis,nel_x,nel_y);
    elseif reduction == 4
        U_basis2 = sum_square_velocities(U_basis,nel_x, nel_y);
    elseif reduction == 5
        U_basis2 = U_basis; 
    elseif reduction == 6
        U_basis1 = generate_symmetric_base(U_basis,nel_x,nel_y);
        ind2_trim = 1:(nel_x/2):(nel_x*nel_y/2);
        U_basis1(:,ind2_trim) = [];
        U_basis2 = U_basis1;
    elseif reduction == 10      % trim out elements in Lithosphere
        elementsOmega1 = InfoMesh.list_Omega1; 
        U_basis2 = U_basis; 
        U_basis2(:,elementsOmega1) = [];
    elseif reduction == 11     % trim out elements in Lithosphere and those in last column
        % define last column 
        elem_last_col_tot = nel_x:nel_x:nel_x*nel_y;
        % retrieve elements in the Lithosphere
        elementsOmega1 = InfoMesh.list_Omega1; 
        % determine those elements in last column not in the Lithosphere
        elem_last_col = elem_last_col_tot(ismember(elem_last_col_tot,elementsOmega1)==0);
        % add them both
        elements2trim = [elementsOmega1; elem_last_col'];
        U_basis2 = U_basis; 
        U_basis2(:,elements2trim) = [];
    elseif reduction == 16 % group elements by four as reduction==4 and trim out last column
        U_basis2 = sum_square_velocities(U_basis,nel_x, nel_y);
        elem2trim = (nel_x/2):(nel_x/2):size(U_basis2,2);
        U_basis2(:,elem2trim) = [];
    elseif reduction == 17      % this groups by four elements as red.==4 and takes out those completely in Omega1
        elements2erase = InfoMesh.list_Omega1;
        U_basis2 = square_minus_lithosp(U_basis,nel_x,nel_y,elements2erase);
    end 

    [U,S,V] = svd(U_basis2);
    contador_svd = sum(diag(S)/S(1,1) < 1e-8);

    % figure(1); clf; 
    % semilogy(diag(S)/S(1,1),'k-x');
    % grid on; grid minor; set(gca,'FontSize',14);
    % ylabel('\sigma_i / \sigma_{max}','FontSize',14)
    % xlabel(['i=1,...,',num2str(size(U_basis2,2))],'FontSize',12)
    % hold on; 
    % plot([1 size(U_basis2,2)], [1e-8 1e-8],'r-')
    % str1 = ['/Users/marianofernandez/Library/Mobile Documents/com~apple~CloudDocs/00-PhD/07-Reportes/15-Reporte_convergencia_minimizacion/figures/svd_figures/base',num2str(reduction),'2.pdf'];
    % saveas(gcf,str1)

    if contador_svd ~= 0
        u2 = reshape(u_ref,[size(U_basis2(:,1))]);
        n_col_SVD = min(size(S)) - contador_svd;
        U_basis2 = U(:,1:n_col_SVD);
        alpha_min = U_basis2'*u2;
        velo_min = U_basis2*alpha_min;
        S_ast = S(1:n_col_SVD,1:n_col_SVD);
        V_ast = V(:,1:n_col_SVD);
        invS = diag(1./diag(S_ast));         
        delta_min = V_ast * invS * alpha_min;  % it's really alpha_min, but it's not used afterwards 
    else
        alpha_min = (U_basis2'*U_basis2) \ (U_basis2' * u_ref(:));
        velo_min = U_basis2*alpha_min; 
    end
        
norm_dif = norm(velo_min - u_ref(:)) / norm(u_ref(:));
velo_min = reshape(velo_min,size(InfoMesh.X_v));
final_size = size(U_basis2,2);
if do_plots.plot == 1
    dasp1 = (InfoMesh.fin_x - InfoMesh.ini_x) > 2*(InfoMesh.fin_y - InfoMesh.ini_y);
    plots123(velo_min,InfoMesh.X_v,norm_dif, final_size,do_plots,dasp1)
end
    

end 


function plots123(velo_min, Xv, error, num_comp, plots_do,daspect123)
    addpath '01_plots'
    param1.parameters = 0; 
    sc = 1; figure_num = 11; 
    plot_velocities(velo_min, sc,figure_num, plots_do.LABx, plots_do.LABy,Xv,param1)
    if daspect123 == 1
        daspect([1 0.5 1])
    else
        daspect([1 1 1])
    end 
    error1 = round(error,2);
    title(['Least-Squares approximation with ',num2str(num_comp),' components'],'FontSize',13)
    annotation('textbox', [0.13, 0.55, 0.28, 0.07], 'BackgroundColor', 'w', 'String', ['||u^{ref} - u_{LS}|| / || u^{ref} || = ',num2str(error1)], ...
    'FontSize', 12);
end

function U_b_sym = generate_symmetric_base(U_basis,nel_x,nel_y)
% generate_symmetric_base - Reduces the original basis to its half by summing the
% symmetric components with respect to the Y-Y axis.
% INPUTS:
%   U_basis - Matrix of basis velocity components (size: m x n)
%   nel_x - Number of elements in the X direction (must be even)
%   nel_y - Number of elements in the Y direction
%
% OUTPUT:
%   U_b_sym - Symmetric basis matrix

    % Validate that nel_x is even
    if mod(nel_x, 2) ~= 0
        error('nel_x must be an even number for symmetric reduction.');
    end
    U_b_sym = zeros(size(U_basis,1), size(U_basis,2)/2);
    % Remember U_basis is constructed from top to bottom in Y-Direction

    for ii = 1:nel_y
        ini_left = (ii-1)*nel_x+1;
        fin_left = ini_left + nel_x/2 - 1;
        left_side = ini_left:1:fin_left;
        ini_right = fin_left+1;
        fin_right = ii*nel_x; 
        right_side = flip(ini_right:1:fin_right);
        U_new = U_basis(:,left_side) + U_basis(:,right_side);
        ind2save = (ii-1)*(nel_x/2)+1:ii*(nel_x/2);
        U_b_sym(:,ind2save) = U_new; 
    end 

end

function U_squares = sum_square_velocities(U_basis,nel_x,nel_y)
% sum_square_velocities - Reduces the original basis to its quarter by
% summing two elements in X- and two elements in Y-directions
% INPUTS:
%   U_basis - Matrix of basis velocity components (size: m x n)
%   nel_x - Number of elements in the X direction (must be even)
%   nel_y - Number of elements in the Y direction
% OUTPUT:
%   U_squares - Symmetric basis matrix


    % Validate that nel_x and nel_y are even
    if mod(nel_x, 2) ~= 0 || mod(nel_y, 2) ~= 0
        error('nel_x and nel_y must be even numbers for square reduction.');
    end
    
    % Initialize the reduced basis matrix
    U_squares = zeros(size(U_basis, 1), (size(U_basis, 2) / 4));

    % Counter to keep track of the column index in the reduced matrix
    counter = 0;

    % Loop over half the elements in the Y direction
    for jj = 1:(nel_y / 2)
        % Loop over half the elements in the X direction
        for ii = 1:(nel_x / 2)
            % Calculate indices for the four components in the original matrix
            base_index = (jj - 1) * 2 * nel_x + (ii - 1) * 2 + 1;
            base_ind1 = base_index;
            base_ind2 = base_ind1 + 1;
            vert_ind1 = base_ind1 + nel_x;
            vert_ind2 = base_ind2 + nel_x;
          %  [base_ind1 base_ind2 vert_ind1 vert_ind2]

            % Sum the four components
            U_new = U_basis(:, vert_ind1) + U_basis(:, vert_ind2) + ...
                    U_basis(:, base_ind1) + U_basis(:, base_ind2);

            % Increment the counter and save the new summed component
            counter = counter + 1;
            U_squares(:, counter) = U_new;
        end
    end
end

function U_squares = square_minus_lithosp(U_basis,nel_x,nel_y,lithosph_elem)
% sum_square_velocities - Reduces the original basis to its quarter by
% summing two elements in X- and two elements in Y-directions and keeps
% only those who are not in list_elem
% INPUTS:
%   U_basis - Matrix of basis velocity components (size: m x n)
%   nel_x - Number of elements in the X direction (must be even)
%   nel_y - Number of elements in the Y direction
% OUTPUT:
%   U_squares - Symmetric basis matrix


    % Validate that nel_x and nel_y are even
    if mod(nel_x, 2) ~= 0 || mod(nel_y, 2) ~= 0
        error('nel_x and nel_y must be even numbers for square reduction.');
    end
    
    % Initialize the reduced basis matrix
    U_squares = zeros(size(U_basis, 1), (size(U_basis, 2) / 4));

    % Counter to keep track of the column index in the reduced matrix
    counter = 0;

    % flag for trimming out elements not needed. 
    flag1 = zeros((size(U_basis,2)/4),1);

    % Loop over half the elements in the Y direction
    for jj = 1:(nel_y/2)
        % Loop over half the elements in the X direction
        for ii = 1:(nel_x/2)
            % Calculate indices for the four components in the original matrix
            base_index = (jj - 1) * 2 * nel_x + (ii - 1) * 2 + 1;
            base_ind1 = base_index;
            base_ind2 = base_ind1 + 1;
            vert_ind1 = base_ind1 + nel_x;
            vert_ind2 = base_ind2 + nel_x;
            elements_considered = [base_ind1 base_ind2 vert_ind1 vert_ind2];
            how_many_in_lithosphere = sum(ismember(elements_considered,lithosph_elem));
            
            % Sum the four components
            U_new = U_basis(:, vert_ind1) + U_basis(:, vert_ind2) + ...
                    U_basis(:, base_ind1) + U_basis(:, base_ind2);

            % Increment the counter and save the new summed component
            counter = counter + 1;
            flag1(counter) = how_many_in_lithosphere == 4;        % all the elements belong to Omega1
            U_squares(:, counter) = U_new;
        end
    end
    %last_col = nel_x/2:nel_x/2:size(U_squares,2);   % last column of components (in Omega1, Omega2, and Gamma)
    %flag1(last_col) = ones(length(last_col),1);
    U_squares(:,logical(flag1)) = [];
end