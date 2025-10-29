function LS_val = LevelSet(X_point,InfoLAB,varargin)
% OUTPUT : 
% LS_val: value of the level set
% INPUT: 
%  X_point : point to evaluate the value of the level set
% InfoLAB: structure where the interface points is
% varargin: Plotting purposes only 

% The interface
LABx=InfoLAB.LABx / (InfoLAB.maxDepth/1000); 
LABy = ((InfoLAB.maxDepth/1000) + InfoLAB.LABy) / (InfoLAB.maxDepth/1000);

% Build polygons
% top boundary
x1 = [-1 2];
y1 = 100;
% lower boundary == interface                                                
x_int = LABx';
y_int = LABy';
% polygons
px = [x1(1) x_int x1(2) x1(2) x1(1) x1(1)]';
py = [y_int(1) y_int  y_int(end) y1 y1 y_int(1)]';

% Calculate distance
[~,Z] = dsearchn([x_int' y_int'],X_point);
in = inpolygon(X_point(:,1),X_point(:,2),px,py); % define those points inside the polygon 
Z(in) = -Z(in);     % set to minus the value 
LS_val = -Z;        % set all to minus the value 

% %% Plot
plot_LS = varargin{1}.plotLS;
if plot_LS == 1
    nel_x = varargin{1}.nel_x; 
    nel_y = varargin{1}.nel_y; 
    if size(X_point,1) == (2*nel_x+1)*(2*nel_y+1) 
        n_grid_x = 2*nel_x+1;
        n_grid_y = 2*nel_y + 1;
    elseif size(X_point,1) == (nel_x+1)*(nel_y+1)
        n_grid_x = nel_x+1;
        n_grid_y = nel_y+1;
    else
        error('REVISAR')
    end
    
    X11 = reshape(X_point(:,1),[n_grid_x n_grid_y]);
    Y11 = reshape(X_point(:,2),[n_grid_x n_grid_y]);
    figure(1); clf; 
    tolerance_plot = varargin{1}.tolPlot;
    Contours = [-1:0.2:-0.2 -0.1 -tolerance_plot:0.001:tolerance_plot 0.1 0.15 0.3];
    %Contours = [-0.1 -tolerance_plot:0.001:tolerance_plot 0.1];
    [C,h] = contour(X11,Y11,reshape(LS_val,size(X11)),Contours);
    set(gca,'FontSize',12)
    h.LineWidth = 1.2;
    clabel(C,h,Contours,'FontSize',12);
    ne_x = varargin{1}.nel_x; 
    ne_y = varargin{1}.nel_y;
    xticks(0:1/ne_x:1)
    yticks(0:1/ne_y:1)
    grid on
    hold on
    plot(LABx,LABy,'b--','LineWidth',2)
    axis tight equal
    title('Level set function','FontSize',14)
end 

end 