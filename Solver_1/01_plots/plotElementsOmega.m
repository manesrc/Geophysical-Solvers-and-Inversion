function plotElementsOmega(Omega,InfoMesh)
    X = InfoMesh.X;
    T = InfoMesh.T;
    figure(120)
    hold on
    % elements
    n_elem_tot = length(Omega);
    for iCol = 1:n_elem_tot
      ce = Omega(iCol);
     Te = T(ce,:);
     Xe = X(Te,:);%*InfoProblem.L_ref;
     addElem(Xe,{'-'});
    end

end


function addElem(X,arg)
   %ix = [1 2 3 4 1 5 6 7 8 5 6 2 3 7 8 4]; % malla lineal
   ix = [1 2 3 4 1];
   plot(X(ix,1),X(ix,2),arg{:})
end