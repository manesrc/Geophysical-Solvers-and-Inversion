function T_2_nd2 = compute_all_dof_inf(vel_01,K_fem2,N2,G2_matrix,f_fem2,vn2,matG1,matG2,T_1_nd,DOF1,DOF2,DOF_Q,PosMat,plot_flux,plot_temp,InfoLAB,InfoProblem,InfoMesh)


% create v_0 and B_0
% New Mat 2
MAT2_new2 = K_fem2+N2+vel_01*G2_matrix;

if InfoMesh.nne == 9
    eleType1D = 0; % line
    nne_1D = 3; % quadratic
elseif InfoMesh.nne == 4
    eleType1D = 0; % line
    nne_1D = 2; % quadratic
end 
[xgp,wgp] = quadrature(eleType1D,nne_1D);
[N,~,~] = shapeFunctions(eleType1D,nne_1D,xgp);
x_length_elem = (1/InfoMesh.nel_x);

% to consider the g_mean integrate over DOF_Q [DOF in model bottom]
line_int_vect = zeros(length(DOF_Q),1);
nel_x = InfoMesh.nel_x;
for jj = 1:nel_x
    if nne_1D == 3
        Te = [2*jj-1 2*jj 2*jj+1];
    else
        Te = [jj jj+1];
    end 
    g_e = zeros(nne_1D,1);
    for kk = 1:size(wgp,2)
        g_e = g_e + N(kk,:)'*wgp(kk) * (x_length_elem/2);
    end
    line_int_vect(Te,1) = line_int_vect(Te,1) + g_e;
end 
g_mean = InfoProblem.q2*line_int_vect;
g_mean_long = [g_mean; zeros(length(DOF2)-length(g_mean),1)];
% vect for v_0 with DOF2 length
vect2_2 = f_fem2 + vn2;


method = 1;
relationship = 1e-3;

fluxOmega1 = matG1*T_1_nd;
fluxZeroOmega2 = matG2 * (MAT2_new2(DOF2,DOF2)\vect2_2(DOF2));
fluxOmega2_gmean = matG2 * (MAT2_new2(DOF2,DOF2)\g_mean_long);
% v_0
v_0 = -fluxOmega1 - fluxZeroOmega2 - fluxOmega2_gmean;


%
pos_flux_eval = PosMat*InfoMesh.X;
matrix_svd = [];

for ii = 1:length(DOF_Q)
    % D_mat enables certain DOF to be different than zero or not
    D_mat2 = zeros(length(DOF_Q),1);
    D_mat2(ii) = 1;
    D_mat = [D_mat2; zeros(length(DOF2)-size(D_mat2,1),size(D_mat2,2))];
    % B_0: 
    B_0 = matG2 * (MAT2_new2(DOF2,DOF2)\D_mat);

    area_cont = g_mean(ii)/InfoProblem.q2;
    
    [g_perturbation,~] = obtain_g2(B_0,v_0,area_cont,method,relationship);

    g_2nn = 1;  % debería ser g_perturbation
    g_2_new = D_mat*g_2nn;

    vect2_new = f_fem2(DOF2) + g_2_new + vn2(DOF2) + g_mean_long;
    
    % CALCULATE T_2
    t2 = tic;
    T_2_nd2 = MAT2_new2(DOF2,DOF2)\vect2_new;
    toc(t2)
    
    % 2D Plot 
    if plot_temp == 1
        g_2_new2plot = g_mean;
        x2plot_qs = linspace(0,1,length(g_2_new2plot));
        num_plot = 113; 
        varargin.perturbation = D_mat2*g_2nn;
        plotBothSub_adim(T_1_nd,T_2_nd2,InfoMesh,InfoLAB,InfoProblem,DOF1,DOF2,g_2_new2plot,x2plot_qs,num_plot,varargin);
        saveas(gcf,['99_figures\2d_temp_plot_DOF',num2str(ii),'.pdf'])
    end 

    % plot jump of flux
    %% NODAL JUMP OF FLUX VALUES
    if plot_flux == 1
        vect2_mean = f_fem2(DOF2) + g_mean_long + vn2(DOF2);
        T_2_mean = MAT2_new2(DOF2,DOF2)\vect2_mean;
        fluxOmega2_gmean = matG2*T_2_mean;
        fluxOmega2 = matG2*T_2_nd2;
        matrix_svd = [matrix_svd fluxOmega2-fluxOmega2_gmean];
        flux_jump = fluxOmega1 + fluxOmega2;  
        figure(55); clf;
        scatter(pos_flux_eval(:,1),flux_jump,'bx');  
        fontSize55 = 16;
        set(gca,'FontSize',fontSize55)
        grid on
        hold on 
        scatter(pos_flux_eval(:,1),fluxOmega1,50,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',0.4);%'ro'); 
        scatter(pos_flux_eval(:,1),-fluxOmega2,50,'^','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0],'MarkerFaceAlpha',0.4);%,'g^');
        scatter(pos_flux_eval(:,1),-fluxOmega2_gmean,25,'v','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[220 220 220]/660,'MarkerFaceAlpha',0.4);%,'g^');
        fluxes123 = [fluxOmega1; -fluxOmega2; flux_jump];  
        
        plot([0 1],[0 0],'k-')
        axis([0 1 round(min(fluxes123),2) round(max(fluxes123),2)]);
        xlabel('Points in \Gamma','FontSize',fontSize55); 
        ylabel('\Delta flux = f_1 + f_2','FontSize',fontSize55)
        legend('(k_1 \nabla T_1-k_2 \nabla T_2)\cdotn_1','k_1\nablaT_1\cdotn_1','-k_2\nablaT_2\cdotn_2','-k_2\nablaT_2(mean)\cdotn_2','Location','southoutside','NumColumns',3,'FontSize',fontSize55)
        saveas(gcf,['99_figures\jump_flux_DOF',num2str(ii),'.pdf'])
    else
        matrix_svd = [];
    end 
end 


 % D_mat enables certain DOF to be different than zero or not
D_mat2 = eye(size(DOF_Q,1));
D_mat = [D_mat2; zeros(length(DOF2)-size(D_mat2,1),size(D_mat2,2))];
% B_0: 
B_0 = matG2 * (MAT2_new2(DOF2,DOF2)\D_mat);

area_cont = diag(g_mean);

[g_perturbation,~] = obtain_g2(B_0,v_0,area_cont,method,relationship);

g_2_new = D_mat*g_perturbation;

vect2_new = f_fem2(DOF2) + g_2_new + vn2(DOF2) + g_mean_long;

% CALCULATE T_2
t2 = tic;
T_2_nd2 = MAT2_new2(DOF2,DOF2)\vect2_new;
toc(t2)

ii = 99;

g_2_new2plot = g_mean;
x2plot_qs = linspace(0,1,length(g_2_new2plot));
num_plot = 113; 
varargin.perturbation = D_mat2*g_perturbation;
plotBothSub_adim2c(1523*T_1_nd,1523*T_2_nd2,InfoMesh,InfoLAB,InfoProblem,DOF1,DOF2,g_2_new2plot,x2plot_qs,num_plot,varargin);
saveas(gcf,['99_figures\2d_temp_plot_DOF',num2str(ii),'.pdf'])

vect2_mean = f_fem2(DOF2) + g_mean_long + vn2(DOF2);
T_2_mean = MAT2_new2(DOF2,DOF2)\vect2_mean;
fluxOmega2_gmean = matG2*T_2_mean;
fluxOmega2 = matG2*T_2_nd2;
flux_jump = fluxOmega1 + fluxOmega2;  

figure(55); clf;
pos_flux_eval = 660*pos_flux_eval;
scatter(pos_flux_eval(:,1),flux_jump,'bx');  
fontSize55 = 18;
set(gca,'FontSize',fontSize55)
grid on
hold on 
scatter(pos_flux_eval(:,1),fluxOmega1,50,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',0.4);%'ro'); 
scatter(pos_flux_eval(:,1),-fluxOmega2,50,'^','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0],'MarkerFaceAlpha',0.4);%,'g^');
scatter(pos_flux_eval(:,1),-fluxOmega2_gmean,25,'v','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[220 220 220]/660,'MarkerFaceAlpha',0.4);%,'g^');
fluxes123 = [fluxOmega1; -fluxOmega2; flux_jump; -fluxOmega2_gmean];  

plot([min(pos_flux_eval(:,1)) max(pos_flux_eval(:,1))],[0 0],'k-')
axis([min(pos_flux_eval(:,1)) max(pos_flux_eval(:,1)) round(min(fluxes123),2) round(max(fluxes123),2)]);
xlabel('Points in \Gamma','FontSize',fontSize55); 
ylabel('flux','FontSize',fontSize55)
legend('(k_1 \nabla T_1-k_2 \nabla T_2)\cdotn_1','k_1\nablaT_1\cdotn_1','-k_2\nablaT_2\cdotn_2','-k_2\nablaT_2(mean)\cdotn_2','Location','southoutside','NumColumns',3,'FontSize',fontSize55)
saveas(gcf,['99_figures\jump_flux_DOF',num2str(ii),'.pdf'])

figure(56); clf;
scatter(pos_flux_eval(:,1),flux_jump,'bx');  
fontSize55 = 18;
set(gca,'FontSize',fontSize55)
grid on
hold on 
scatter(pos_flux_eval(:,1),fluxOmega1,50,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',0.4);%'ro'); 
scatter(pos_flux_eval(:,1),-fluxOmega2,50,'^','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0],'MarkerFaceAlpha',0.4);%,'g^');
%scatter(pos_flux_eval(:,1),-fluxOmega2_gmean,25,'v','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[220 220 220]/660,'MarkerFaceAlpha',0.4);%,'g^');
fluxes123 = [fluxOmega1; -fluxOmega2; flux_jump];  
plot([0 660],[0 0],'k-')
axis([0 660 round(min(fluxes123),2) round(max(fluxes123),2)]);
xlabel('Points in \Gamma','FontSize',fontSize55); 
ylabel('flux','FontSize',fontSize55)
legend('flux_1+flux_2','flux_1','-flux_2','Location','southoutside','NumColumns',3,'FontSize',fontSize55)
text(330,0.1,'Performance = 10.9 %','FontSize',fontSize55)
saveas(gcf,['99_figures\jump_flux_DOF',num2str(ii+1),'.pdf'])


if isempty(matrix_svd) ~= 1
    [~,D,~] = svd(matrix_svd);
    vect2plot = diag(D)/D(1,1);
    figure(85); clf;
    semilogy(vect2plot,'k-x')
    set(gca,'FontSize',14)
    xlabel('Eigenmodes','FontSize',14)
    ylabel('Eigenvalues / 1st. eigval. ( \lambda_i / \lambda_{1,1} )','FontSize',14)
    axis tight
    grid on 
    grid minor
    saveas(gcf,['99_figures\SVD.pdf'])
end 





end