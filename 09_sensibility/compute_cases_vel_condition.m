function T_2_nd2 = compute_cases_vel_condition(vel_01,position,K_fem2,N2,G2_matrix,f_fem2,vn2,matG1,matG2,T_1_nd,DOF1,DOF2,DOF_Q,PosMat,InfoLAB,InfoProblem,InfoMesh)


% create v_0 and B_0
% New Mat 2
MAT2_new2 = K_fem2+N2+vel_01*G2_matrix;
% to consider the g_mean integrate over DOF_Q [DOF in model bottom]
line_int_vect = zeros(length(DOF_Q),1);
nel_x = InfoMesh.nel_x;
[xgp,wgp] = quadrature(0,3);
[N,~,~] = shapeFunctions(0,3,xgp);
x_length_elem = (1/InfoMesh.nel_x);
for jj = 1:nel_x
    Te = [2*jj-1 2*jj 2*jj+1];
    g_e = zeros(3,1);
    for kk = 1:size(wgp,2)
        g_e = g_e + N(kk,:)'*wgp(kk) * (x_length_elem/2);
    end
    line_int_vect(Te,1) = line_int_vect(Te,1) + g_e;
end 
g_mean = InfoProblem.q2*line_int_vect;
g_mean_long = [g_mean; zeros(length(DOF2)-length(g_mean),1)];
% vect for v_0 with DOF2 length
vect2_2 = f_fem2 + vn2;



% D_mat enables certain DOF to be different than zero or not
D_mat2 = zeros(length(DOF_Q),1);
if position == 1
    D_mat2(InfoMesh.nel_x+1) = 1;
elseif position == 2
    D_mat2(7) = 1;
elseif position == 3
    D_mat2(end-6) = 1;
elseif position == 4
    D_mat2 = zeros(length(DOF_Q),3);
    D_mat2(7,1) = 1;
    D_mat2(InfoMesh.nel_x+1,2) = 1;
    D_mat2(end-6,3) = 1;
end 


D_mat = [D_mat2; zeros(length(DOF2)-size(D_mat2,1),size(D_mat2,2))];


fluxOmega1 = matG1*T_1_nd;
fluxZeroOmega2 = matG2 * (MAT2_new2(DOF2,DOF2)\vect2_2(DOF2));
fluxOmega2_gmean = matG2 * (MAT2_new2(DOF2,DOF2)\g_mean_long);
v_0 = -fluxOmega1 - fluxZeroOmega2 - fluxOmega2_gmean;
B_0 = matG2 * (MAT2_new2(DOF2,DOF2)\D_mat);




method = 2;
relationship = 1e-3;

area_cont = g_mean(find(D_mat2)~=0)/InfoProblem.q2;

[g_perturbation,q_perturbation] = obtain_g2(B_0,v_0,area_cont,method,relationship);

g_2_new = D_mat*g_perturbation;

vect2_new = f_fem2(DOF2) + g_2_new + vn2(DOF2) + g_mean_long;

t2 = tic;
T_2_nd2 = MAT2_new2(DOF2,DOF2)\vect2_new;
toc(t2)


% 2D Plot 
g_2_new2plot = g_mean;
x2plot_qs = linspace(0,1,length(g_2_new2plot));
num_plot = 113; 
varargin.perturbation = D_mat2*g_perturbation;


plotBothSub_adim(T_1_nd,T_2_nd2,InfoMesh,InfoLAB,InfoProblem,DOF1,DOF2,g_2_new2plot,x2plot_qs,num_plot,varargin);
%disp(['Solution g = ',num2str(g_solution2)])


%% NODAL JUMP OF FLUX VALUES
fluxOmega1 = matG1*T_1_nd;
fluxOmega2 = matG2*T_2_nd2;
flux_jump = fluxOmega1 + fluxOmega2;

pos_flux_eval = PosMat*InfoMesh.X;

figure(55); clf;
scatter(pos_flux_eval(:,1),flux_jump,'bx');  
fontSize55 = 16;
set(gca,'FontSize',fontSize55)
grid on
hold on 
scatter(pos_flux_eval(:,1),fluxOmega1,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',0.4);%'ro'); 
scatter(pos_flux_eval(:,1),-fluxOmega2,'^','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0],'MarkerFaceAlpha',0.4);%,'g^');
fluxes123 = [fluxOmega1; -fluxOmega2; flux_jump];
%fluxes123 = [flux_jump; jump_flux_Dolbow];

if position==1
    text1 = 'central DOF ≠ 0';
elseif position == 2
    text1 = 'left DOF ≠ 0';
elseif position == 3
    text1 = 'right DOF ≠ 0';
elseif position == 4
    text1 = 'all together';
end 

%text(0.05,0.195,['Min. elem. ratio to check flux cont. = ',num2str(min_ratio)],'FontSize',14)
text(0.25,0.18,text1,'FontSize',fontSize55)
plot([0 1],[0 0],'k-')
axis([0 1 round(min(fluxes123),2) round(max(fluxes123),2)]);
xlabel('Points in \Gamma','FontSize',fontSize55); 
ylabel('\Delta flux = f_1 + f_2','FontSize',fontSize55)
legend('(k_1 \nabla T_1-k_2 \nabla T_2)\cdotn_1','k_1\nablaT_1\cdotn_1','-k_2\nablaT_2\cdotn_2','Location','southoutside','NumColumns',3,'FontSize',fontSize55)
%legend('residue = (k_1 - k_2) \nabla T n_1','k_1 \nabla T n_1','-k_2 \nabla T n_2','Location','southoutside', ...
  %  'NumColumns',3,'FontSize',11)

norm_jump = flux_jump'*flux_jump;
normflux1 = fluxOmega1'*fluxOmega1;
residue1 = norm_jump/normflux1;
perf = 100*round((1-residue1),3);
disp(['The performance of the code (Delta flux = flux_1 + flux_2) is ',num2str(perf),'% true'])





end