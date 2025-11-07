function LS_val = LevelSet(X_point,InfoLAB,varargin)
% this function evaluates the level set function for the LAB interface
% OUTPUT : 
% LS_val: value of the level set
% INPUT: 
%  X_point : point to evaluate the value of the level set
% InfoLAB: structure where the interface points is
% varargin: Plotting purposes only 

% The interface
npoints2addLAB = 1000;
[LABx,LABy] = increase_data_LAB(InfoLAB,npoints2addLAB);


% Build polygons
% top boundary
x1 = [min(LABx)-max(LABx) 2*max(LABx)];
y1 = 2*max(X_point(:,2));
% lower boundary == interface                                                
if min(LABx) == LABx(1)
    x_int = LABx';
    y_int = LABy';
else
    x_int = flip(LABx)';
    y_int = flip(LABy)';
end 
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
    Y2 = Y11 - max(X_point(:,2));
    figure(1); clf; 
    [C,h] = contour(X11,abs(Y2),reshape(LS_val,size(X11)));
    font_size = 12; 
    set(gca,'YDir','reverse','FontSize',font_size)
    h.LineWidth = 1.2;
    clabel(C,h,'FontSize',font_size);
    grid on
    hold on
    LABy2 = max(X_point(:,2)) - LABy; 
    plot(LABx,LABy2,'b--','LineWidth',2)
    xlabel('X-Axis [m]','FontSize',font_size)
    ylabel('Depth [m]','FontSize',font_size)
    axis tight equal
    title('Level set function','FontSize',font_size+2)
end 

end 