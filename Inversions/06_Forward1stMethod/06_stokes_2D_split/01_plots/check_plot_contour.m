function check_plot_contour(fig_num,magnitude,X)
x_d = length(unique(X(:,1)));
y_d = length(unique(X(:,2)));
[xr, yr] = meshgrid(unique(X(:,1)), unique(X(:,2)));
figure(fig_num); clf;
T_r = reshape(magnitude, [x_d y_d] );
[C,h] = contourf(xr,yr,T_r');
%colormap jet
colorbar
h.LineWidth = 1.2;
clabel(C,h,'FontSize',12);
c2 = colorbar;
c2.Limits = [min(magnitude) max(magnitude)]; 

end
