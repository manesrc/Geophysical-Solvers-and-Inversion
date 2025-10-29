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
addpath '06_stokes_2D/00_utils'
addpath '06_stokes_2D/01_plots'
% Mesh for pressure (create one)
elemType = InfoMesh.elemType;
nenPre = 4; % create a lineal interpolation 
if nenPre == InfoMesh.nne
    X = InfoMesh.X;
    T = InfoMesh.T;
    XP = X;
    TP = T;
    % generar malla para velocidad nne = 9
    InfoMesh.NewMesh = 1;  % en
    InfoMesh.nne_new = 9;
    [X,T] = CreaMalla_rectangulo_MTF(InfoMesh);
    % cambiar las velocidades a la nueva malla
    [Temp12,Temp22,LS_mesh] = convTempMesh(X,T,TP,Temp1,Temp2,InfoMesh);
else 
    X = InfoMesh.X;
    T = InfoMesh.T;
    InfoMesh.NewMesh = 1;  % en
    InfoMesh.nne_new = 4;
    [XP,TP] = CreaMalla_rectangulo_MTF(InfoMesh);
    Temp12 = Temp1;
    Temp22 = Temp2;
    LS_mesh = InfoMesh.LS_mesh;
end



%% material properties

% viscosity 
mu_ref = InfoMaterial.mu_ref;
% density [top_domain; bottom_domain]
rho_ref = InfoMaterial.rho_ref ;
% gravity
InfoMaterial.L_ref = InfoProblem.maxDepth;
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
    LABx = plot_up.LAB(:,1)/660;
    LABy = ((L_ref/1000) + plot_up.LAB(:,2))/660;
    % plot velocities
    figure(4); clf
    subplot(121)
    hold on
%    scale_velo = velo;
    quiver(X(:,1),X(:,2),velo2(:,1),velo2(:,2),1,'k');   % *(L_ref/1000)
    plot(LABx,LABy,'r--')
    box on
    axis equal tight
    title('velocities')

    
    % plot pressures
    subplot(122)
    tri = delaunay(XP(:,1),XP(:,2));
    pres = full(pres);
    trisurf(tri,XP(:,1)*(L_ref/1000),XP(:,2)*(L_ref/1000),pres*(1e-6));
    grid on
    axis tight
    title('pressures')
    xlabel('X-Axis')
    ylabel('Y-Axis')
    zlabel('Pressure in [MPa]')


    figure(120);
    hold on
    quiver(X(:,1),X(:,2),velo2(:,1),velo2(:,2),1,'k');   % *(L_ref/1000)

    if plot_up.parameters == 1
        LABData = plot_up.LAB;
        %LABData.LABy = plot_up.LAB(:,2);
        Contours = [0.19:0.2:1 1.0:0.01:1.5]; % these are curves ploted
        title_fig = 'Temperature [K] (T) ';
        makecontour_check2(temp,xgp22,101,InfoMesh.nel_x,InfoMesh.nel_y,LABData,Contours,title_fig)
        Contours = [0.8:0.05:1.3]; % these are curves ploted
        title_fig = '\rho (\rho/3300 kg/m^3)';
        makecontour_check2(rho/rho_ref,xgp22,102,InfoMesh.nel_x,InfoMesh.nel_y,LABData,Contours,title_fig)
        Contours = [16:0.5:25];
        title_fig = 'log_{10} (\mu) ';
        makecontour_check2(log10(mu),xgp22,103,InfoMesh.nel_x,InfoMesh.nel_y,LABData,Contours,title_fig)
    end 
    
    LABx = LABx*660;
    LABy = (LABy-1)*660;

    % plot the quiver with less arrows
    sc = 5;
    arguments1 = 1:sc:size(X,1);
    
    sf1 = figure(88); clf;
    quiver(X(arguments1,1)*660,660*(X(arguments1,2))-660,velo2(arguments1,1),velo2(arguments1,2),1,'k')
    axis tight
    set(gca,'FontSize',18);
    %h1.YTick = 0:-100:-600;
    sf1.CurrentAxes.YTickLabel = {'600','500','400','300','200','100','0'};
    xlabel('X [km]','FontSize',18)
    ylabel('Depth [km]','FontSize',18)
    daspect([1 1 1])
    hold on
    plot(LABx,LABy,'r-','LineWidth',3)
    str1 = '$$ \Gamma_{_{\textsf{LAB}}} $$';
    text(480,-200,str1,'Color','red','FontSize',16,'Interpreter','latex') 
    %legend('velocities','LAB','FontSize',18)
    max_abs_u = 365*100*24*3600*max(abs(velo));
    str2 = ['max(|u|) = [',num2str(round(max_abs_u(1),2)),' ',num2str(round(max_abs_u(2),2)),'] cm/yr'];
    ht = text(50,-600,str2) ;
    set(ht, 'color','k','backgroundcolor','w','EdgeColor','k','FontSize',12)
end



end