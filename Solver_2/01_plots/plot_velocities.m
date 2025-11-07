function plot_velocities(velo,sc,num_fig,PropChangeX,PropChangeY,X,param)
    % function to plot velocities as quiver plot
    % inputs: 
    % velocity field (velo): field to plot,
    % scale (sc): scale for quiver plot (less arrows for increasing sc)
    % figure number (num_fig): figure number for the plot
    % proper change in X (PropChangeX): x-coordinates of the proper change line (LAB)
    % proper change in Y (PropChangeY): y-coordinates of the proper change line (LAB)
    % coordinates (X): mesh coordinates
    % param: structure with additional parameters (optional title, parameters to plot as densities and viscosities)
    
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
    font_size = 20; 
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
    plot(PropChangeX,PropChangeY,'r-','LineWidth',3)
    %legend('','velocities','\Delta proper.','FontSize',font_size ,'location','best')


    % Check if 'title' exists in param
    if isfield(param, 'title') && ~isempty(param.title)
        title(param.title, 'FontSize', font_size + 2)
    end


    if param.parameters == 1
        vect1 = [max(PropChangeX), max(PropChangeY), max(X(:,1)), max(X(:,2))];
        L_refe1 = max(vect1);
        mu_used = param.mu_ref * param.mu_used;
        rho_used = param.rho_ref * param.rho_used;
        xgp22 = param.xgp22;
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