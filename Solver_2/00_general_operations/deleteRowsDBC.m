function [K, f] = deleteRowsDBC(K, f, DBCmatrix)
% deleteRowsDBC: Modifies the global stiffness matrix K and force vector f to apply 
% Dirichlet BC
% outputs:
%   K - modified global stiffness matrix    
%   f - modified global force vector
% inputs:
%   K - global stiffness matrix
%   f - global force vector
%   DBCmatrix - matrix where the first column contains the row indices to be modified

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


    