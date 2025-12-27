function check_plot_contour(fig_num,magnitude,X)
% check_plot_contour: Plots a filled contour of the provided magnitude
% over the mesh defined by coordinates X.
% INPUT:
%   fig_num   : Figure number for the plot
%   magnitude : Vector of magnitudes to be contoured at each node
%   X         : Nx2 matrix of nodal coordinates (X(:,1) = x
%               coordinates, X(:,2) = y coordinates)

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
