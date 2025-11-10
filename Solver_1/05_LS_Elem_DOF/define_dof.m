function [DOF1,DOF2,DOF_inter,DOF_Q] = define_dof(InfoMesh)
% define the degrees of freedom active in each sub-problem
% INPUT:
%   InfoMesh: structure containing mesh information
% OUTPUT:
%   DOF1: vector with the DOF active in subdomain Omega1
%   DOF2: vector with the DOF active in subdomain Omega2
%   DOF_inter: vector with the DOF active in the interphase
%   DOF_Q: vector with the DOF in the bottom boundary (to apply flux BC)

T = InfoMesh.T;                              % connectivity matrix
list1 = InfoMesh.list1;                     % elements in Omega1
list2 = InfoMesh.list2;                     % elements in Omega2
list_cut = InfoMesh.list_cut;      % elements well-crossed by interface
list_edge1 = InfoMesh.list_edge1;
list_edge2 = InfoMesh.list_edge2;

if isempty(list_cut) ~= 1 && isempty(list_edge1) ~= 1
    list1_complete = [list1; list_cut(:,1); list_edge1(:,1)];
    list2_complete = [list2; list_cut(:,1); list_edge2(:,1)];
    intersected_DOF_arista = DOF_in_arista(T,list_edge1,list_edge2);
    DOF_wellCross = T(list_cut(:,1),:);
elseif isempty(list_cut) ~= 1
    list1_complete = [list1; list_cut(:,1)];
    list2_complete = [list2; list_cut(:,1)];
    intersected_DOF_arista = [];
    DOF_wellCross = T(list_cut(:,1),:);
elseif isempty(list_edge1) ~= 1
    list1_complete = [list1; list_edge1(:,1)];
    list2_complete = [list2; list_edge2(:,1)];
    intersected_DOF_arista = DOF_in_arista(T,list_edge1,list_edge2);
    DOF_wellCross = [];
else
    error('REVISAR DOF')
end 
DOF1_repeated = T(list1_complete,:);
DOF2_repeated = T(list2_complete,:);

DOF1 = unique(DOF1_repeated(:));
DOF2 = unique(DOF2_repeated(:));
DOF_inter = unique([intersected_DOF_arista; DOF_wellCross(:)]);

X = InfoMesh.X;
DOF_Q = find(X(:,2) == min(X(:,2)));

end


% this function retrieves the DOF in the interphase and considers only those DOF in 
% list_edge coinciding with the interphase, disregarding the others
function DOF_illcrossed2retain = DOF_in_arista(T,list_edge1,list_edge2)
numnodes_in_arista = sqrt(size(T,2));
elem_illcrossed1 = list_edge1(:,1); 
elem_illcrossed2 = list_edge2(:,1);
num_elem_illCrossed = size(elem_illcrossed1,1);
DOF_illcrossed2retain1 = zeros(num_elem_illCrossed,numnodes_in_arista);
for qqm = 1:num_elem_illCrossed
    contador1 = 0;
    elem_Omega1 = elem_illcrossed1(qqm);
    elem_Omega2 = elem_illcrossed2(qqm);
    Te1 = T(elem_Omega1,:);
    Te2 = T(elem_Omega2,:);
    for pp = 1:length(Te1)
        if ismember(Te1(pp),Te2)
            contador1 = contador1+1;
            DOF_illcrossed2retain1(qqm,contador1) = Te1(pp);
        end 
    end 
end 
DOF_illcrossed2retain = unique(DOF_illcrossed2retain1(:));
end 