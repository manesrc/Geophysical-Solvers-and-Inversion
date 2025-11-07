function iso_Temp = plot_results(TempObj,do_plot,results,fig_num, InfoMesh,InfoLAB,varargin)
    % function to plot temperature results
    % inputs:
    % TempObj: temperature value to plot the isotherm
    % do_plot: flag to show (1) or not (0) the plot (useful for iterative procedures)
    % results: temperature field to plot (in Kelvin)
    % fig_num: figure number for the plot 
    % InfoMesh: structure with mesh information (X-coordinates)
    % InfoLAB: structure with LAB information (LABx, LABy)
    % varargin: optional structure with additional parameters (e.g., title)
    
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
    contours = [20 200:300:1100 1300:100:1700];
    [C,h] = contourf(xr,abs(y2),T_r',contours); 
    set(gca,'YDir','reverse')
    font_size = 16;
    set(gca,'FontSize',font_size)
    xlabel('X-Axis [m]','FontSize',font_size)
    ylabel('Depth [m]','FontSize',font_size)
    hold on


    [C1,~] = contour(xr,abs(y2),T_r',[1 1] * TempObj_inC,'g--','LineWidth',3);
    % set color limitst and colormap
    colormap(jet);
    caxis([20 1600]);
    iso_Temp = unique(C1','rows')';
    find_iso12 = iso_Temp(1,:) == TempObj_inC;
    iso_Temp(:,find_iso12) = [];
    h.LineWidth = 1.2;
    clabel(C,h,'FontSize',font_size-2);
    c2 = colorbar;
    c2.Limits = [20 1600]; 
    ylabel(c2,'Temperatures [ºC]','FontSize',font_size)
    LABy2 = max(X(:,2)) - InfoLAB.LABy;
    plot(InfoLAB.LABx,LABy2,'b--','LineWidth',2)


    if isempty(varargin) == 0
        str = varargin{1}.title;    
        title(str,'FontSize',font_size+2)
    end 

    if do_plot == 0
        set(0,'DefaultFigureVisible','on')
    end 

    daspect([1 1 1])


end 