function [velo,pres] = ComputeStokesProblem(Temp1,Temp2,InfoMesh,InfoProblem,InfoMaterial,plot_up)
% the idea of this function is compute the velocity of the mantle flow and
% the pressure from Stokes Problem
% OUTPUTS: 
%                   u: velocity as a vector of DOF in mesh x 2 directions [nDOF_X,2]
%                   p: pressure as a vector of DOF in mesh x 1 [nDOF_X,1]
% INPUTS:
%                   T1: temperature distribution in Omega_1
%                   T2: temperature distribution in Omega_2
%                   LS: The level set values to decide whether the point of interest belongs to 
%                           Omega_1 or Omega_2
%                   InfoMesh: characteristics of the mesh, the mesh, etc
%                   InfoProblem: characteristics of the problem, lengths 
%                   InfoMaterial: description to obtain material characteristics such to compute 
%                                       rho and mu, etc

%% Meshes
addpath '06_Forward1stMethod/06_stokes_2D_split/00_utils'
addpath '06_Forward1stMethod/06_stokes_2D_split/01_plots'
% Mesh for pressure (create one)
elemType = InfoMesh.elemType;
X = InfoMesh.X_v / InfoProblem.L_ref;   % has to be dimensionless in Stokes
T = InfoMesh.T_v;
XP = InfoMesh.X / InfoProblem.L_ref; % has to be dimensionless in Stokes
TP =  InfoMesh.T; 
[Temp12,Temp22,LS_mesh] = convTempMesh(X,T,TP,Temp1,Temp2,InfoMesh);

%% material properties

% viscosity 
mu_ref = InfoMaterial.mu_ref;
% density [top_domain; bottom_domain]
rho_ref = InfoMaterial.rho_ref ;
% gravity
InfoMaterial.L_ref = InfoProblem.L_ref;
InfoMaterial.gravity = InfoMaterial.gravity_units*( (InfoMaterial.L_ref^3) * (rho_ref^2)/ (mu_ref^2) );

%% Build matrices
[K,G,f,mu,temp,rho,xgp22] = makeGlobalMatrices(elemType,X,T,XP,TP,InfoMaterial,Temp12,Temp22,LS_mesh,InfoMesh,InfoProblem,plot_up);

%% Boundary conditions
ksize = size(K,1);
[Adir,bdir] = bcFreeSlip(X,ksize,XP);

% build the global system
gsize = size(G,1);
Z01 = sparse(gsize,gsize);
Z0cc = sparse( size(Adir,1),size(Adir,1) );

Atot = [[K G'; G Z01; Adir] [Adir'; Z0cc] ];
btot = [f; zeros(gsize,1); bdir];

% solve it!
aux = Atot\btot;

% back to units:
mu_ref = InfoMaterial.mu_ref;
L_ref = InfoMaterial.L_ref;
rho_ref = InfoMaterial.rho_ref;

% extract the velocity 
velo1 = aux(1:ksize) * (mu_ref/(L_ref*rho_ref)); % -> back to [m/s]
velo2 = reshape(velo1,2,[])';

% extract pressures
%pres1 = aux(ksize+1:ksize+gsize) * ;
pres = aux(ksize+1:ksize+gsize) * ( (mu_ref^2) / ( (L_ref^2)*rho_ref )  );  % -> back to [Pa]

% convert the results into a linear mesh
if InfoMesh.nne == 4
    velo = convertQuad2Lin_velo(velo2,XP,TP,T);
else
    velo = velo2;
end 


%% Post process
if plot_up.fig == 1
    LABx = plot_up.LAB(1,:);
    LABy = plot_up.LAB(2,:);
    % % plot velocities
    % figure(4); clf
    % subplot(121)
    % hold on
    % quiver(InfoProblem.L_ref * X(:,1),InfoProblem.L_ref * X(:,2),velo2(:,1),velo2(:,2),1,'k');   % *(L_ref/1000)
    % plot(LABx,LABy,'r--')
    % box on
    % axis equal tight
    % title('velocities')
    % % plot pressures
    % subplot(122)
    % tri = delaunay(XP(:,1),XP(:,2));
    % pres = full(pres);
    % trisurf(tri,XP(:,1)*(L_ref/1000),XP(:,2)*(L_ref/1000),pres*(1e-6));
    % grid on
    % axis tight
    % title('pressures')
    % xlabel('X-Axis')
    % ylabel('Y-Axis')
    % zlabel('Pressure in [MPa]')
    % figure(120);
    % hold on
    % quiver(InfoProblem.L_ref *X(:,1),InfoProblem.L_ref *X(:,2),velo2(:,1),velo2(:,2),1,'k');   % *(L_ref/1000)
   
    % plot the quiver with less arrows
    sc = 5;
    %arguments1 = 1:sc:size(X,1);
    
    sf1 = figure(88); 
    clf;
    % plot the quiver with less arrows
    arguments1 = 1:sc:size(X,1);
    % to plot in depth: 
    newY = InfoProblem.L_ref * (max(X(:,2)) - X(arguments1,2));
    veloY = -velo2(arguments1,2);
    quiver(InfoProblem.L_ref * X(arguments1,1),newY,velo2(arguments1,1),veloY,1,'k','LineWidth',1);
    set(gca,'YDir','reverse','FontSize',18 )
    axis tight
    %xlabel('X-Axis [m]','FontSize',font_size)
    %ylabel('Depth [m]','FontSize',font_size)
    %daspect([1 1 1])
    %quiver(InfoProblem.L_ref *X(arguments1,1),InfoProblem.L_ref*(max(X(arguments1,2))-X(arguments1,2)),velo2(arguments1,1),velo2(arguments1,2),1,'k')
    %axis tight
    %set(gca,'YDir','reverse','FontSize',18);
    xlabel('X [km]','FontSize',18)
    ylabel('Depth [km]','FontSize',18)
    daspect([1 0.5 1])
    hold on
    plot(InfoProblem.L_ref * LABx,InfoProblem.L_ref * LABy,'r-','LineWidth',3)
    str1 = '$$ \Gamma_{_{\textsf{LAB}}} $$';
    text(InfoProblem.L_ref/2,0.5*InfoProblem.L_ref,str1,'Color','red','FontSize',16,'Interpreter','latex') 
    %legend('velocities','LAB','FontSize',18)
    max_abs_u = 365*100*24*3600*max(abs(velo));
    str2 = ['max(|u|) = [',num2str(round(max_abs_u(1),2)),' ',num2str(round(max_abs_u(2),2)),'] cm/yr'];
    ht = text(InfoProblem.L_ref*0.2,InfoProblem.L_ref*0.2, str2) ;
    set(ht, 'color','k','backgroundcolor','w','EdgeColor','k','FontSize',12)
end

if plot_up.parameters == 1
    L_refe1 = InfoProblem.L_ref; 
    rho_used = rho; mu_used = mu; font_size = 12; 
    PropChangeX = plot_up.LAB(1,:);
    PropChangeY = plot_up.LAB(2,:);
    figure(131); clf;
    scatter3(L_refe1*xgp22(:,1), L_refe1 * abs(max(xgp22(:,2)) - xgp22(:,2) ) , rho_used);
    set(gca,'YDir','reverse')
    hold on; 
    plot3(L_refe1 *PropChangeX,L_refe1 *PropChangeY,max(rho_used)*ones(size(PropChangeY)),'r--','LineWidth',3)
    axis tight
    xlabel('X-Axis [m]','FontSize',font_size)
    ylabel('Depth [m]','FontSize',font_size)
    zlabel('\rho [kg/m^3]','FontSize',font_size)
    figure(132); clf;
    scatter3(L_refe1*xgp22(:,1), L_refe1 * abs(max(xgp22(:,2)) - xgp22(:,2) ) , log10(mu_used));
    hold on
    plot3(L_refe1 * PropChangeX,L_refe1 *PropChangeY,max(log10(mu_used))*ones(size(PropChangeY)),'r--','LineWidth',3)
    %plot3(PropChangeX,PropChangeY,max(log10(mu_used))*ones(size(PropChangeY,1),1),'r--','LineWidth',3)
    axis tight
    xlabel('X-Axis [m]','FontSize',font_size)
    ylabel('Depth [m]','FontSize',font_size)
    zlabel('log(\mu) [log(Pa s)]','FontSize',font_size)
end 


end