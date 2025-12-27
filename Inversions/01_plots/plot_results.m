function iso_Temp = plot_results(TempObj,do_plot,results,fig_num, InfoMesh,InfoLAB,varargin)
% plot_results: Plots the temperature results over the mesh along with
% the LAB position.
% INPUT:
%   TempObj   : Temperature value defining the LAB (in Kelvin)
%   do_plot   : Flag to control plot visibility (1 = visible, 0 = invisible)
%   results   : Vector of temperature results at each node (in Kelvin)
%   fig_num   : Figure number for the plot
%   InfoMesh  : Structure containing mesh information including:
%               - X: Node coordinates
%   InfoLAB   : Structure containing LAB position information including:
%               - LABx: X-coordinates of the LAB line
%               - LABy: Y-coordinate of the LAB line
%   varargin  : Optional structure with additional plot settings:
%               - title: Title for the plot
%

X = InfoMesh.X;
x_d = length(unique(X(:,1)));
y_d = length(unique(X(:,2)));
results_inC = results - 273; 
TempObj_inC = TempObj - 273; 
T_r = reshape(results_inC, [x_d y_d] );
[xr, yr] = meshgrid(unique(X(:,1)), unique(X(:,2)));
y2 = yr - max(X(:,2));

if do_plot == 0
    set(0,'DefaultFigureVisible','off')
end 
h1 = figure(fig_num); clf;
contours = [0 200:300:1100 TempObj_inC 1300:20:ceil(max(results_inC))];
[C,h] = contourf(xr,abs(y2),T_r',contours); 
set(gca,'YDir','reverse')
font_size = 12;
set(gca,'FontSize',font_size)
xlabel('X-Axis [m]','FontSize',font_size)
ylabel('Depth [m]','FontSize',font_size)
hold on

[C1,~] = contour(xr,abs(y2),T_r',[1 1] * TempObj_inC,'g--','LineWidth',3);
iso_Temp = unique(C1','rows')';
find_iso12 = iso_Temp(1,:) == TempObj_inC;
iso_Temp(:,find_iso12) = [];
colormap jet
h.LineWidth = 1.2;
clabel(C,h,'FontSize',10);
c2 = colorbar;
c2.Limits = [300 1500]; 
ylabel(c2,'Temperatures [ºC]','FontSize',font_size)
LABy2 = max(X(:,2)) - InfoLAB.LABy;
plot(InfoLAB.LABx,LABy2,'b-','LineWidth',2)


if isempty(varargin) == 0
    str = varargin{1}.title;    
    title(str,'FontSize',font_size+2)
end 

if do_plot == 0
    set(0,'DefaultFigureVisible','on')
end 

daspect([1 1 1])


end 