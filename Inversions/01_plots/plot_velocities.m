function plot_velocities(velo,sc,num_fig,PropChangeX,PropChangeY,X,param)
% this function plots the velocity field
% INPUT:
%   velo          : Nx2 matrix of nodal velocities (velo(:,1) = vx
%                   components, velo(:,2) = vy components)
%   sc            : scaling factor for quiver plot (to reduce number of arrows)
%   num_fig       : Figure number for the plot
%   PropChangeX   : X-coordinates of the property change line
%   PropChangeY   : Y-coordinates of the property change line
%   X             : Nx2 matrix of nodal coordinates (X(:,1) = x
%                   coordinates, X(:,2) = y coordinates)
%   param         : structure containing parameter information including:
%                   - parameters: flag indicating if parameters are to be plotted
%                   - mu_ref: reference viscosity
%                   - mu_used: viscosity scaling factors at Gauss points
%                   - rho_ref: reference density
%                   - rho_used: density scaling factors at Gauss points
%                   - xgp22: Gauss point coordinates


    velo1 = full(velo); % [m/s]
x1 = unique(X(:,1)); 
y1 = unique(X(:,2));
z1 = sqrt(velo1(:,1).^2+velo1(:,2).^2); 
z2 = reshape(z1,[length(x1) length(y1)])*100*365*24*60*60; % [cm/yr]
X1 = reshape(X(:,1),[length(x1) length(y1)]);
Y1 = reshape(X(:,2),[length(x1) length(y1)]);
Y2 = Y1 - max(X(:,2));

PropChangeY = max(X(:,2)) - PropChangeY;

% plot the quiver with less arrows
figure(num_fig); clf;
% plot surface of the vector magnitude
[c1,h1] = contourf(X1,abs(Y2),z2);
set(gca,'YDir','reverse');
font_size = 12; 
shading interp
h1.LineStyle = 'none';
colormap('summer')
c2 = colorbar;
ylabel(c2,'||velocities||_{L_2} [cm/yr]','FontSize',font_size )
hold on
% plot the quiver with less arrows
arguments1 = 1:sc:size(X,1);
% to plot in depth: 
newY = max(X(:,2)) - X(arguments1,2);
veloY = -velo(arguments1,2);
quiver(X(arguments1,1),newY,velo(arguments1,1),veloY,1,'k','LineWidth',1);
set(gca,'FontSize',font_size )
axis tight
xlabel('X-Axis [m]','FontSize',font_size)
ylabel('Depth [m]','FontSize',font_size)
daspect([1 1 1])
hold on
plot(PropChangeX,PropChangeY,'b--','LineWidth',3)
legend('','velocities','\Delta proper.','FontSize',font_size ,'location','best')

if param.parameters == 1
    vect1 = [max(PropChangeX), max(PropChangeY), max(X(:,1)), max(X(:,2))];
    L_refe1 = max(vect1);
    mu_used = param.mu_ref * param.mu_used;
    rho_used = param.rho_ref * param.rho_used;
    xgp22 = param.xgp22;
    % % create a contourplot for parameters
    % % Step 1: Sort the data based on x and y coordinates
    % [sortedX, sortedIndex] = sort(xgp22(:,1));
    % sortedY = xgp22(sortedIndex,2);
    % sortedParameterValue_rho = rho_used(sortedIndex);
    % sortedParameterValue_mu = mu_used(sortedIndex);
    % % Step 2: Create a structured grid using meshgrid
    % [Xgrid, Ygrid] = meshgrid(unique(sortedX), unique(sortedY));
    % % Step 3: Interpolate parameter values at the grid points
    % interpolatedParameterValue_mu = griddata(sortedX, sortedY, sortedParameterValue_mu, Xgrid, Ygrid);
    % interpolatedParameterValue_rho = griddata(sortedX, sortedY, sortedParameterValue_rho, Xgrid, Ygrid);
    % Xgrid_units = L_refe1* Xgrid;
    % Ygrid_units = L_refe1* Ygrid;
    % Ygrid_units_depth = Ygrid_units - max(X(:,2));
    % % Step 4: Plot the contour
    % figure(152); clf; 
    % contours_mu = [1e19; 1e20; 1e21; 1e22; 1e23; 1e24];
    % contourf(Xgrid_units, abs(Ygrid_units_depth), interpolatedParameterValue_mu,contours_mu);
    % set(gca,'YDir','reverse')
    % daspect([1 1 1])
    % colorbar;
    % hold on; 
    % plot(PropChangeX,PropChangeY,'r--','LineWidth',3)
    % set(gca,'FontSize',font_size)
    % title('\mu','FontSize',font_size)
    % xlabel('X-Axis [m]','FontSize',font_size)
    % ylabel('Depth [m]','FontSize',font_size)
    % title('\mu at Gauss Stokes points (plot w/ Matlab interp.)','FontSize',font_size)
    % figure(151); clf; 
    % contours_rho  = 3000:50:5000;
    % [C,h] = contourf(Xgrid_units, abs(Ygrid_units_depth), interpolatedParameterValue_rho, contours_rho);
    % set(gca,'YDir','reverse')
    % clabel(C,h,'FontSize',font_size-2);
    % daspect([1 1 1])
    % hold on; 
    % plot(PropChangeX,PropChangeY,'r--','LineWidth',3)
    % colorbar;
    % set(gca,'FontSize',font_size)
    % title('\rho at Gauss Stokes points (plot w/ Matlab interp.)','FontSize',font_size)
    % xlabel('X-Axis [m]','FontSize',font_size)
    % ylabel('Depth [m]','FontSize',font_size)
    % set(gca,'FontSize',12)
    figure(131); clf;
    scatter3(L_refe1*xgp22(:,1), L_refe1 * abs(max(xgp22(:,2)) - xgp22(:,2) ) , rho_used);
    set(gca,'YDir','reverse')
    hold on; 
    plot3(PropChangeX,PropChangeY,max(rho_used)*ones(size(PropChangeY,1),1),'r--','LineWidth',3)
    axis tight
    xlabel('X-Axis [m]','FontSize',font_size)
    ylabel('Depth [m]','FontSize',font_size)
    zlabel('\rho [kg/m^3]','FontSize',font_size)
    figure(132); clf;
    scatter3(L_refe1*xgp22(:,1), L_refe1 * abs(max(xgp22(:,2)) - xgp22(:,2) ) , log10(mu_used));
    hold on
    plot3(PropChangeX,PropChangeY,max(log10(mu_used))*ones(size(PropChangeY,1),1),'r--','LineWidth',3)
    axis tight
    xlabel('X-Axis [m]','FontSize',font_size)
    ylabel('Depth [m]','FontSize',font_size)
    zlabel('log(\mu) [log(Pa s)]','FontSize',font_size)
end 

end 