function plot_execute(InfoLAB,L_ref,execution_plot)
% plot showing the LAB distribution if execution_plot == 0 doesn't plot it 
% input:
% InfoLAB : structure containing the LAB data
% L_ref   : reference length for axis scaling
% execution_plot : flag to execute the plot (1: plot, 0: no plot)
LABx = InfoLAB.LABx;
LABy = L_ref/1000+InfoLAB.LABy;
%L_ref = InfoLAB.maxDepth/1000;
    if execution_plot ~= 0
        figure(249); clf;
        plot(LABx,LABy,'b-*','MarkerSize',1); 
        axis([0 L_ref 0 L_ref])
        grid on ; 
        grid minor; 
        xlabel('x-axis','FontSize',12); 
        ylabel('y-axis','FontSize',12)
        title('isotherm defining the LAB','FontSize',12)
    end 
end 