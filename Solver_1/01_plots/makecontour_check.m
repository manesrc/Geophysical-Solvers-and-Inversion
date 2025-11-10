function makecontour_check(var,pointsXY,plot_num,ne_x,ne_y,LABData,Contours,title_fig)
% receives a vectorial variable and a set of ponts and plots a contour plot
% INPUTS :
%   var         : variable to plot (vector)
%   pointsXY    : coordinates of the points where var is defined (matrix with 2 columns)
%   plot_num    : figure number to use
%   ne_x        : number of elements in the X direction
%   ne_y        : number of elements in the Y direction
%   LABData     : LAB curve data (matrix with 2 columns)
%   Contours    : contour levels to plot (vector)
%   title_fig   : title of the figure (string)

x_points = pointsXY(:,1);
y_points = pointsXY(:,2);
num_points = length(var);
points_grid = sqrt(num_points);

LABx = LABData(:,1)/660;
LABy = (660 + LABData(:,2))/660;


matrix_plot = zeros(points_grid);
matrix_coordX = zeros(points_grid);
matrix_coordY = zeros(points_grid);

ngp = num_points/(ne_x*ne_y) ;
ngp_dir = sqrt(ngp);


cont_matCol = 1;
cont_matRow = 1;


for ii = 1:ne_y
    for jj = 1:ne_x
        cont_ini = (ii-1) * (num_points/ne_x) +  (jj-1) * ngp + 1;
        cont_fin = (ii-1) * (num_points/ne_x) + jj*ngp;
        points_elementX = reshape(x_points(cont_ini:cont_fin),[ngp_dir ngp_dir])';
        points_elementY = reshape(y_points(cont_ini:cont_fin),[ngp_dir ngp_dir])';
        var_data = reshape(var(cont_ini:cont_fin),[ngp_dir ngp_dir])';
        contCol = cont_matCol:(cont_matCol+ngp_dir-1);
        contRow = cont_matRow:(cont_matRow+ngp_dir-1);
        matrix_plot(contRow,contCol) = var_data;
        matrix_coordX(contRow,contCol) = points_elementX;
        matrix_coordY(contRow,contCol) = points_elementY;
        cont_matCol = cont_matCol + ngp_dir;
    end
    cont_matRow = cont_matRow +ngp_dir;
    cont_matCol = 1;
end

%Contours = 1505*[0.19:0.2:8 0.99 1.0:0.01:1.5]; % these are curves ploted

figure(plot_num);
clf;
hold on
[Cs,h] = contourf(matrix_coordX,matrix_coordY,matrix_plot,Contours);
colormap('jet')
plot(LABx,LABy,'b--','LineWidth',1.5)
set(gca,'FontSize',12)
h.LineWidth = 1.2;
clabel(Cs,h,Contours,'FontSize',14)
title(title_fig,'FontSize',14)


figure(plot_num+100);
clf;
hold on
scatter3(x_points,y_points,var,5,var)
xlabel('X-Axis [-]','FontSize',12)
ylabel('Y-Axis [-]','FontSize',12)
zlabel(title_fig,'FontSize',12)
colormap jet
colorbar
plot3(LABx,LABy,min(var(:))*ones(size(LABx)),'r--','LineWidth',1.5)
set(gca,'FontSize',12)
h.LineWidth = 1.2;
clabel(Cs,h,Contours,'FontSize',14)
title(title_fig,'FontSize',14)
%contour(matrix_coordX,)