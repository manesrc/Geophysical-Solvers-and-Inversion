function [K, f] = deleteRowsDBC(K, f, DBCmatrix)
% DELETEROWSDBC Modifies the stiffness matrix K and the load vector f
% to apply Dirichlet boundary conditions specified in DBCmatrix.
% INPUT:
%   K          - Stiffness matrix
%   f          - Load vector
%   DBCmatrix  - Matrix where each row contains [rowIndex, prescribedValue]
% OUTPUT:
%   K          - Modified stiffness matrix
%   f          - Modified load vector
% Dirichlet boundary conditions -------------------------------------------
rowsToDelete = DBCmatrix(:,1);
nOfRowsToDelete = length(rowsToDelete);
DBCmatrix = DBCmatrix(:,1:2);

for kRow = 1:nOfRowsToDelete
    iRow = rowsToDelete(kRow);
    f = f - DBCmatrix(kRow,2)*K(:,iRow);
end

K(rowsToDelete,:) = 0;
K(:,rowsToDelete) = 0;
for kRow = 1:nOfRowsToDelete
    iRow = rowsToDelete(kRow);
    K(iRow,iRow) = 1;
end
f(rowsToDelete) = DBCmatrix(:,2);
end 


    