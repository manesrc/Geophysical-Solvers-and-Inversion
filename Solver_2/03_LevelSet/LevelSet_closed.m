function LS_val = LevelSet_closed(X_point,polygon,plot_LS)
% this function evaluates the level set function for closed polygons (used to compute velocity basis)
% OUTPUT : 
% LS_val: value of the level set
% INPUT: 
%  X_point : point to evaluate the value of the level set
% polygon: vertices of the area where the different properties are applied
% plot_LS: Plotting purposes only 


Xdense = addpoints2Polygon(polygon,50);
% Closed polygons
x_int = Xdense(:,1);
y_int = Xdense(:,2);



% Calculate distance
[~,Z] = dsearchn([x_int y_int],X_point);
in = inpolygon(X_point(:,1),X_point(:,2),x_int,y_int); % define those points inside the polygon 
Z(in) = -Z(in);     % set to minus the value 
LS_val = -Z;        % set all to minus the value 

% %% Plot
if plot_LS.plotLS == 1
    nel_x = plot_LS.nel_x; 
    nel_y = plot_LS.nel_y; 
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
    tolerance_plot = plot_LS.tolPlot;
    L_ref1 = [1 max(max(X_point(:,1), max(X_point(:,2))))];
    fact = [max(L_ref1) == 1 max(L_ref1)~= 1]* [1 max(max(X_point(:,1), max(X_point(:,2))))]';
    Contours = fact * [-1:0.2:-0.2 -0.1 -tolerance_plot:0.001:tolerance_plot 0.1 0.15 0.3];
    %Contours = [-0.1 -tolerance_plot:0.001:tolerance_plot 0.1];
    [C,h] = contour(X11,Y11,reshape(LS_val,size(X11)),Contours);
    set(gca,'FontSize',12)
    h.LineWidth = 1.2;
    clabel(C,h,Contours,'FontSize',12);
    xticks(min(X_point(:,1)):(max(X_point(:,1))-min(X_point(:,1)))/nel_x:max(X_point(:,1)))
    yticks(min(X_point(:,2)):(max(X_point(:,2))-min(X_point(:,2)))/nel_y:max(X_point(:,2)))
    grid on
    hold on
    plot(x_int,y_int,'b--','LineWidth',2)
    axis tight equal
    title('Level set function','FontSize',14)
end 

end 