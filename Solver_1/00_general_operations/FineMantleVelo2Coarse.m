function [u_coarse,p_coarse] = FineMantleVelo2Coarse(results_fine,InfoMesh_fine,InfoMesh_coarse)
% this function receives a mantle flow velocity field of a fine mesh and adapts the results to a 
% coarse mesh to avoid these variations modify results

u_coarse = zeros(size(InfoMesh_coarse.X));
p_coarse = zeros(size(InfoMesh_coarse.X,1),1);
u_fine = results_fine.u_mantle;
p_fine = results_fine.p_mantle;
X_coarse = InfoMesh_coarse.X;
posit_mesh_mat = zeros(size(X_coarse,1),1);
X_fine = InfoMesh_fine.X;

for ii = 1:size(u_coarse,1)
    refer_point = X_coarse(ii,:);
    [~,posit_in_fine_mesh] = ismember(refer_point,X_fine, 'rows');
    posit_mesh_mat(ii) = posit_in_fine_mesh;
    u_coarse(ii,:) = u_fine(posit_in_fine_mesh,:);
    p_coarse(ii) = p_fine(posit_in_fine_mesh);
end 




%% plot to check 

% % plot velocities
% figure(4), clf
% subplot(121)
% hold on
% quiver(X_fine(:,1),X_fine(:,2),u_fine(:,1),u_fine(:,2),1,'k');   % *(L_ref/1000)
% box on
% axis equal tight
% title('original velocity')
% subplot(122)
% quiver(X_coarse(:,1),X_coarse(:,2),u_coarse(:,1),u_coarse(:,2),1,'k');   % *(L_ref/1000)
% axis equal tight
% title('adapted velocity')
% 
% figure(120);
% hold on
% quiver(X_fine(:,1),X_fine(:,2),u_fine(:,1),u_fine(:,2),1,'k');   % *(L_ref/1000)
% 
% figure(121);
% hold on
% quiver(X_coarse(:,1),X_coarse(:,2),u_coarse(:,1),u_coarse(:,2),1,'k');   % *(L_ref/1000)


end 