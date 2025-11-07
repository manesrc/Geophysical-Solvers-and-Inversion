function check_plot_contour(fig_num,magnitude,X,InfoLAB,title2)
% check_plot_contour: function to plot contour plot of magnitude data
% inputs:   
%   fig_num     : figure number
%   magnitude   : magnitude data to plot    
%   X           : coordinates of mesh points
%   InfoLAB     : structure with LAB information (optional)
%   title2      : title of the plot (optional)

    x_d = length(unique(X(:,1)));
    y_d = length(unique(X(:,2)));
    [xr, yr] = meshgrid(unique(X(:,1)), unique(X(:,2)));
    figure(fig_num); clf;
    maxY = max(X(:,2));
    T_r = reshape(magnitude, [x_d y_d] );
    [C,h] = contourf(xr,maxY-yr,T_r');
    set(gca,'YDir','reverse')
    %colormap jet
    colorbar
    h.LineWidth = 1.2;
    clabel(C,h,'FontSize',12);
    c2 = colorbar;
    if min(magnitude) ~= max(magnitude)
        c2.Limits = [min(magnitude) max(magnitude)]; 
    end 


    if isempty(InfoLAB) == 0
        hold on 
        plot(InfoLAB.LABx,maxY-InfoLAB.LABy,'r--','LineWidth',2)
    end 

    if isempty(title2) == 0
        title(title2,'FontSize',14)
    end 

end
