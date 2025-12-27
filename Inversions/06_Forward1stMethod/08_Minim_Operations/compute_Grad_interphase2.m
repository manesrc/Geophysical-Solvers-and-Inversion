function [PosMatrix,PosMatrix1,matrix_Gradiente1,PosMatrix2,PosMatrix21,matrix_Gradiente2] = compute_Grad_interphase2(DOF1,DOF2,figurePointsInterface,InfoProblem,InfoMesh)

% Mesh related information
T = InfoMesh.T;
X = InfoMesh.X; %/InfoProblem.L_ref;
elemType = InfoMesh.elemType;
nne = InfoMesh.nne;
nel_x = InfoMesh.nel_x; 
LS_mesh = InfoMesh.LS_new;

% find which elements crossed by the interface have one element above the
% lithosphere and one below the asthenosphere
list_cut = InfoMesh.list_cut;
logi1= ones(size(list_cut,1),1);
elem_Omega1 = [];
elem_Omega2 = [];
for ii = 1:size(list_cut,1)
    elem_ii = InfoMesh.list_cut(ii,1);
    Te = T(elem_ii,:);
    Xe = X(Te,:);
    normal1_ii = InfoMesh.list_cut(ii,end-1:end);       % normal to Omega_1 [always by construction]
    if abs(normal1_ii(1)) < abs(normal1_ii(2))        % normal_y > normal_x 
        elem_potentially_Omega1 = elem_ii - (sign(normal1_ii(2) ) ) *  nel_x; % if normal(2) < 1 sums an element line
        elem_potentially_Omega2 = elem_ii + (sign(normal1_ii(2) ) ) * nel_x;
    else            % normal_x > normal_y
        elem_potentially_Omega1 = elem_ii - (sign(normal1_ii(1) ) ) * 1; % if normal(1) < 1 sums one element
        elem_potentially_Omega2 = elem_ii + (sign(normal1_ii(1) ) ) * 1;
    end 
    if (ismember(elem_potentially_Omega1,list_cut(:,1)) == 0) && (ismember(elem_potentially_Omega2,list_cut(:,1)) == 0)
        % if it enters here means there may be elements to compute the gradient 
            logi1(ii) = 0;
            LS_potentially_Omega1 = LS_mesh(T(elem_potentially_Omega1,:));
            LS_potentially_Omega2 = LS_mesh(T(elem_potentially_Omega2,:));
            if (sum(sign(LS_potentially_Omega1)) == 4) && (sum(sign(LS_potentially_Omega2)) == -4)
                info2save1 = [elem_potentially_Omega1 normal1_ii];
                info2save2 = [elem_potentially_Omega2 -normal1_ii];
                elem_Omega1 = [elem_Omega1; info2save1];
                elem_Omega2 = [elem_Omega2; info2save2];
            end 
    end 
end 

% for ii = 1:size(elem_Omega1,1)
%     figure(120); hold on;
%     scatter(X(T(elem_Omega1(ii),:),1),X(T(elem_Omega1(ii),:),2),'bo')
%     scatter(X(T(elem_Omega2(ii),:),1),X(T(elem_Omega2(ii),:),2),'rx')
% end
if isempty(InfoMesh.list_edge1) ~= 1
    list_edge1_mod = [InfoMesh.list_edge1(:,1) InfoMesh.list_edge1(:,end-1:end)]; 
    list_edge2_mod = [InfoMesh.list_edge2(:,1) -InfoMesh.list_edge2(:,end-1:end)]; 
else
    list_edge1_mod = [];
    list_edge2_mod = [];
end

elements2computeGrad1 = [elem_Omega1; list_edge1_mod];
elements2computeGrad2 = [elem_Omega2; list_edge2_mod];

% Obtain S matrix, the one with the value of the shape functions evaluated on points of the 
% the interface where the jump of flux is measured
tot_dof = size(InfoMesh.X,1);
[PosMatrix,PosMatrix1,matrix_Gradiente1] = computeGradientMatrix(elements2computeGrad1,DOF1,InfoProblem.k1,InfoMesh.elemType,T,X);
[PosMatrix2,PosMatrix21,matrix_Gradiente2] = computeGradientMatrix(elements2computeGrad2,DOF2,InfoProblem.k2,InfoMesh.elemType,T,X);
end 

function [PosMatrix,PosMatrix2,G_matrix] = computeGradientMatrix(elements,DOF,conductivity,eleType,T,X)
    % compute the gradient matrix 
    num_elem = size(elements,1);
    nne = size(T,2);
    MatrixGradient = zeros(num_elem,size(X,1)); 
    PosMatrix = zeros(num_elem,size(X,1)); 
    elem_center = [0 0];
    for ii = 1:num_elem
        [N_LABpoints,Nxi_LABp,Neta_LABp] = shapeFunctions(eleType,nne,elem_center);
        % retrieve element connectivities and nodes
        elem_of_interest_jj = elements(ii,1);
        Te = T(elem_of_interest_jj,:);
        Xe = X(Te,:);     
        mat1 = [Nxi_LABp; Neta_LABp];
        Jacob = mat1*Xe; 
        res = Jacob\mat1;
        Nx_qq = res(1,:);
        Ny_qq = res(2,:);
        normal_interest = elements(ii,2:3);
        Contracted_Gradient = [Nx_qq' Ny_qq']*normal_interest'; 
        MatrixGradient(ii,Te) = Contracted_Gradient;  % put the values of the SF in each row
        PosMatrix(ii,Te) = N_LABpoints;  % put the values of the SF in each row
    end 
    G_matrix = conductivity*MatrixGradient(:,DOF);
    PosMatrix2 = PosMatrix(:,DOF); 
end 

% 
% matrix_Position = zeros(nOfPointsFlux,tot_dof);
% % contador 
% matrixRow = 0;
% 
% for jj = 1:nelem2consider           % element by element
%     % intersection points in isop. coord.
%     isop_coord_int1 = intersected_elements(jj,2:3);
%     isop_coord_int2 = intersected_elements(jj,4:5);
%     % normal in the element of interest
%     normal_interest = normal_sign*intersected_elements(jj,end-1:end);
%     % add points between intersection points
%     mid_points_chi = linspace(isop_coord_int1(1),isop_coord_int2(1),nOfPoints_Elem+2);
%     mid_points_eta = linspace(isop_coord_int1(2),isop_coord_int2(2),nOfPoints_Elem+2);
%     iso_coord = [mid_points_chi(2:end-1)' mid_points_eta(2:end-1)'];
%     % obtain shape functions values
%     [N_LABpoints,Nxi_LABp,Neta_LABp] = shapeFunctions(elemType,nne,iso_coord);
% 
%     % retrieve element connectivities and nodes
%     elem_of_interest_jj = intersected_elements(jj,1);
%     Te = T(elem_of_interest_jj,:);
%     Xe = X(Te,:); 
% 
%     % save to the matrix
%     for qq=1:nOfPoints_Elem
%         mat1 = [Nxi_LABp(qq,:); Neta_LABp(qq,:)];
%         Jacob = mat1*Xe; 
%         res = Jacob\mat1;
%         Nx_qq = res(1,:);
%         Ny_qq = res(2,:);
%         Contracted_Gradient = [Nx_qq' Ny_qq']*normal_interest'; 
%         matrix_Gradiente(matrixRow+qq,Te) = Contracted_Gradient;  % put the values of the SF in each row
%         matrix_Position(matrixRow+qq,Te) = N_LABpoints(qq,:);
%     end 
%     matrixRow = matrixRow + nOfPoints_Elem;
%     if figurePointsInterface == 1
%         Xe = InfoMesh.X(Te,:);
%         XPoints = N_LABpoints*Xe;
%         figure(120); hold on; 
%         scatter(XPoints(:,1),XPoints(:,2),20,'kx')
%     end 
% 
% end
% 
% G_matrix = k_diffussivity*matrix_Gradiente(:,DOF_prob);
% matrixS = matrix_Position; %(:,DOF_prob);
% assert(matrixRow == size(G_matrix,1))
% assert(nelem2consider * nOfPoints_Elem == matrixRow)