function [K, f] = deleteRowsDBC(K, f, DBCmatrix)

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


    