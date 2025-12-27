function delta_phi = computeDeltaPhi(dh, Cpoint, ls, phi_e, Xe, InfoMesh)
    % Mesh 
    X = InfoMesh.X;
    % compute derivative of shape functions in isop coord
    center_iso = [0 0];
    [~,Nxi,Neta] = shapeFunctions(InfoMesh.elemType,InfoMesh.nne,center_iso);
    % compute derivative of shape functions in X,Y
    jacob = [Nxi*Xe(:,1)  Nxi*Xe(:,2);
                Neta*Xe(:,1) Neta*Xe(:,2)]; 
    dNxy = jacob\[Nxi; Neta]; 
    dNx_times_phie = dNxy(1,:)*phi_e;
    dNy_times_phie = dNxy(2,:)*phi_e;
    norm_grad = sqrt(dNx_times_phie*dNx_times_phie'+dNy_times_phie*dNy_times_phie');
    dist2Cpoint = sqrt( (X(:,1)-Cpoint(1)).^2 + (X(:,2)-Cpoint(2)).^2 );
    delta_phi = (dh * norm_grad) * exp(-(dist2Cpoint/ls).^2);
end