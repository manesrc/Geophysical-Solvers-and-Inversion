function plot_temp_flux(ii,T2_d,perturbation,fluxOmega2_d,results,InfoMesh,InfoLAB,InfoProblem,varargin)
% post-process temperature and flux plots
% input:
% ii              : figure number
% T2_d           : dimensionless temperature in subdomain 2
% perturbation   : perturbation to add to the gradient for plotting
% fluxOmega2_d  : flux in subdomain 2 (dimensional)
% results        : structure containing results data
% InfoMesh       : structure containing mesh data
% InfoLAB        : structure containing LAB data
% InfoProblem    : structure containing problem data
% varargin       : cell containing additional parameters

%% post-process plots
gradient2plot = results.gradient_T_mean + perturbation;
x2plot_qs = linspace(0,1,length(gradient2plot));
num_plot = ii; 
varargin{1}.perturbation = [];%D_mat2*g_perturbation;
%plotBothSub_adim2c(1523*T_1_nd,1523*T_2_nd,InfoMesh,InfoLAB,InfoProblem,DOF1,DOF2,g_2_new2plot,x2plot_qs,num_plot,varargin);
T_1_d = InfoProblem.T_ref * results.T1nd;
T_2_d = InfoProblem.T_ref * T2_d;

if varargin{1}.restricted == 0 
    plotBothSub_dim_EGU(T_1_d,T_2_d,InfoMesh,InfoLAB,InfoProblem,results.DOF1,results.DOF2,gradient2plot,x2plot_qs,num_plot,varargin);
    saveas(gcf,['99_figures\2d_temp_plot_DOF',num2str(ii),'nelx',num2str(InfoMesh.nel_x),'.pdf'])
    %saveas(gcf,['99_figures\2d_temp_plot_DOF',num2str(ii),'nelx',num2str(InfoMesh.nel_x),'.fig'])
elseif varargin{1}.restricted == 1
    plotBothSub_dim_EGU_rest(T_1_d,T_2_d,InfoMesh,InfoLAB,InfoProblem,results.DOF1,results.DOF2,gradient2plot,x2plot_qs,num_plot,varargin);
    saveas(gcf,['99_figures\2d_temp_plot_DOF_rest',num2str(ii+10),'nelx',num2str(InfoMesh.nel_x),'.pdf'])
    %saveas(gcf,['99_figures\2d_temp_plot_DOF',num2str(ii+10),'nelx',num2str(InfoMesh.nel_x),'.fig'])
end 

fluxOmega1 = results.fluxOmega1_d;
fluxOmega2 = fluxOmega2_d;
fluxOmega2_gmean = results.fluxOmega2_gmean_d;
% dimenisonal
flux_jump = fluxOmega1 + fluxOmega2 ;

figure(ii+15); clf;
pos_flux_eval = results.PosMat*InfoMesh.X;
pos_flux_eval = (InfoProblem.L_ref)*pos_flux_eval;

[pos_flux_orderx,index] = sort(pos_flux_eval(:,1),'ascend');
fluxOmega1_order = fluxOmega1(index);
fluxOmega2_order = fluxOmega2(index);
fluxOmega2gmean_order = fluxOmega2_gmean(index);
jump_flux2 = fluxOmega1_order + fluxOmega2_order;
jump_flux2_gmean = fluxOmega1_order + fluxOmega2gmean_order;

int_flux1 = trapz(pos_flux_orderx,abs(fluxOmega1_order));
int_flux2 = trapz(pos_flux_orderx,abs(fluxOmega2_order));
int_flux2gmean = trapz(pos_flux_orderx,abs(fluxOmega2gmean_order));
int_flux_jump = trapz(pos_flux_orderx,abs(jump_flux2));
int_flux_jump_gmean = trapz(pos_flux_orderx,abs(jump_flux2_gmean));


scatter(pos_flux_eval(:,1),flux_jump,'bx');  
fontSize55 = 16;
set(gca,'FontSize',fontSize55)
grid on
hold on 
scatter(pos_flux_eval(:,1),fluxOmega1,50,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',0.4);%'ro'); 
scatter(pos_flux_eval(:,1),-fluxOmega2,50,'^','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0],'MarkerFaceAlpha',0.4);%,'g^');
scatter(pos_flux_eval(:,1),-fluxOmega2_gmean,25,'v','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[220 220 220]/660,'MarkerFaceAlpha',0.4);%,'g^');
fluxes123 = [fluxOmega1; -fluxOmega2; flux_jump; -fluxOmega2_gmean];  
plot([min(pos_flux_eval(:,1)) max(pos_flux_eval(:,1))],[0 0],'k-')
norm_final_flux = norm(flux_jump);
norm_flux1 = norm(fluxOmega1);
norm_contflux_gmean = norm(fluxOmega1+fluxOmega2_gmean);

perf_minim = 100*round((1 - norm_final_flux/norm_contflux_gmean),3);
ini_flux = round((norm_contflux_gmean/norm_flux1),3);
final_flux = round((norm_final_flux/norm_flux1),3);

str1 = {['Perf. minim. = ',num2str(perf_minim),'%'],['Cont. flux. ini. = ',num2str(ini_flux)],['flux. cont. meas. = ',num2str(final_flux)]};
text(200,min(fluxes123)+1.1*mean(fluxes123),str1,'FontSize',fontSize55)

axis([min(pos_flux_eval(:,1)) max(pos_flux_eval(:,1)) round(min(fluxes123),3) round(max(fluxes123),3)]);
xlabel('Points in \Gamma','FontSize',fontSize55); 
ylabel('flux [W/m^2]','FontSize',fontSize55)
legend('(flux_L-flux_A)','flux_L','-flux_A','-flux_A (mean)','Location','southoutside','NumColumns',2,'FontSize',fontSize55)
saveas(gcf,['99_figures\jump_flux_DOF',num2str(ii),'nelx',num2str(InfoMesh.nel_x),'.pdf'])
%saveas(gcf,['99_figures\jump_flux_DOF',num2str(ii),'nelx',num2str(InfoMesh.nel_x),'.fig'])

end