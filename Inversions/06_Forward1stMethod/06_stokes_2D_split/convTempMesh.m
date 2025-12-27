function [Temp1long,Temp2long,LS_Mesh] = convTempMesh(X,T,TP,Temp1,Temp2,InfoMesh)
% de InfoMesh saco los elementos en Omega 1 y Omega 2. 
% Las temperaturas Temp1 y Temp2 pertenecen a la malla XP, TP -> recorro TP
% uno a uno todos los elementos y calculo temperaturas en los nodos
% restantes
% Temp1long = [nnodos_X,1] <- Temp1 = [nnodos_XP,1]
% Temp2long = [nnodos_X,1] <- Temp2 = [nnodos_XP,1]
nnodosX = size(X,1);
Temp1long = zeros(nnodosX,1); 
Temp2long = zeros(nnodosX,1); 
LS_Mesh = zeros(nnodosX,1); 
% define lists of elements beloging to Omega1 and Omega2 and both

if isempty(InfoMesh.list_edge1) ~= 1
    inOmega1 = [InfoMesh.list1; InfoMesh.list_edge1(:,1)];
    inOmega2 = [InfoMesh.list2; InfoMesh.list_edge2(:,1)];
else
    inOmega1 = InfoMesh.list1;
    inOmega2 = InfoMesh.list2;
end

if isempty(InfoMesh.list_cut) ~= 1
    inBoth = InfoMesh.list_cut(:,1);
else
    inBoth = [];
end 

nElements_TempMesh = size(TP,1);

Xi_gp = [0 -1; 1 0; 0 1; -1 0; 0 0];
[N_lin,~,~] = shapeFunctions(InfoMesh.elemType,InfoMesh.nne,Xi_gp);


for i = 1:nElements_TempMesh
    %DOFs
   DOF_element = TP(i,:);
   DOF_newMesh = T(i,:);
    % Level set values
   LS_ref = InfoMesh.LS_mesh(DOF_element);
   LS_Mesh(DOF_newMesh(1:4)) = LS_ref;
   LS_Mesh(DOF_newMesh(5:end)) = N_lin*LS_ref;
   if ismember(i,inOmega1)
        % Temp_ref + Level set
       Temp_ref = Temp1(DOF_element);       
       Temp1long(DOF_newMesh(1:4)) = Temp_ref;
       Temp1long(DOF_newMesh(5:end)) = N_lin*Temp_ref;
   elseif ismember(i,inOmega2)
       Temp_ref = Temp2(DOF_element);
       Temp2long(DOF_newMesh(1:4)) = Temp_ref;
       Temp2long(DOF_newMesh(5:end)) = N_lin*Temp_ref;
   elseif ismember(i,inBoth)
       % fill T_1
       Temp_ref1 = Temp1(DOF_element);
       Temp1long(DOF_newMesh(1:4)) = Temp_ref1;
       Temp1long(DOF_newMesh(5:end)) = N_lin*Temp_ref1;
       % fill T_2
       Temp_ref2 = Temp2(DOF_element);
       Temp2long(DOF_newMesh(1:4)) = Temp_ref2;
       Temp2long(DOF_newMesh(5:end)) = N_lin*Temp_ref2;
   end 
end 



end