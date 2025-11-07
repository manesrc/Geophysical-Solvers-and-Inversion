function [InfoLAB,InfoMesh,InfoProblem,InfoMaterial,u_mant,p_omega,Temp] = generate_converged_state(cond_dir,InfoLAB,InfoMesh,InfoMaterial,InfoProblem)
% this function gets a LAB and iterates until find a convergence 
% considering the error function and the iterative criterion
% OUTPUTS:
% InfoLAB: stores the converged geometry for the LAB
% InfoMesh: mesh information (with crossed elements by LAB updated)
% InfoProblem: problem information (with updated parameters new values for k in Omega2)
% InfoMaterial: material properties information 
% u_mant: converged velocity field in the mantle
% p_omega: converged pressure field in the whole domain
% Temp: converged temperature field in the whole domain
% INPUTS:
% cond_dir: 1 for Dirichlet condition at the top boundary of Omega
%            2 for Neumann condition at the top boundary of Omega
% InfoLAB: initial geometry for the LAB
% InfoMesh: mesh information
% InfoMaterial: material properties information
% InfoProblem: problem information


warning('this function provides a relevant result for Afonso08 LAB only, the others diverge')
    %% Level-Set for the Gamma_LAB
    addpath '03_LevelSet' 
    % plot-related
    varargin.plotLS = 0; varargin.nel_x = InfoMesh.nel_x; varargin.nel_y = InfoMesh.nel_y;
    InfoMesh.LS_val = LevelSet(InfoMesh.X,InfoLAB,varargin);
    daspect([1 1 1])
    if varargin.plotLS == 1
        saveas(gcf,['99_Figures_Afonso\','1_LS_1','.pdf'])
    end 
    %% Set Stokes properties and run case for Gamma_LAB
    addpath '00_general_operations'
    [T_est,p_est] = estimate_Temp_pres(InfoMesh.X,InfoMesh.nel_x,InfoMesh.nel_y,InfoMesh.LS_val, InfoMaterial,InfoProblem);
    addpath '01_plots'
    varargin.title = 'Temperature estimation';
    plot_results(InfoProblem.T_LAB,1,T_est,5,InfoMesh,InfoLAB,varargin);
    daspect([1 1 1])
    InfoMesh.LS_new = InfoProblem.T_LAB - T_est; % keep Omega1 positive and Omega2 negative
    InfoMesh.LS_val = []; % to ensure its not used anymore
    %% find elements crossed by interphase
    plot1 = 1; 
    tol = 0.001*InfoProblem.T_LAB ; 
    [InfoMesh.list_Omega1,InfoMesh.list_Omega2,InfoMesh.list_cut, InfoMesh.list_edge_elem1, InfoMesh.list_edge_elem2, interphasePoints] = findElemCrossedLS(tol,plot1,InfoMesh, [] ,InfoProblem);
    %% Stokes Problem
    example_plots.fig_velo=1;
    example_plots.parameters = 0; 
    a1 = size(interphasePoints,1) == max(size(interphasePoints));
    maxDepth = max(InfoMesh.X(:,2));
    interphasePoints(2,:) = maxDepth - interphasePoints(2,:);
    if a1 == 0
        example_plots.LAB = [interphasePoints(1,:)'/InfoProblem.L_ref interphasePoints(2,:)'/InfoProblem.L_ref ];    % has to be a two column matrix
    else
        example_plots.LAB = [interphasePoints(1,:)/InfoProblem.L_ref  interphasePoints(2,:)/InfoProblem.L_ref ];  % has to be a two column matrix
    end 
    addpath '05_Stokes_2D'
    [u_mant,p_omega] = ComputeStokesProblem(T_est,p_est, InfoMesh,InfoProblem,InfoMaterial,example_plots);
    daspect([1 1 1])
    if example_plots.fig_velo == 1
        saveas(gcf,['99_Figures_Afonso\2_u_conv1','.pdf'])
    end 
    %% Calculate the temperature field using u_mant
    addpath '04_TempInMesh'
    Temp = solveTemperatureProblem(u_mant,T_est,p_est,cond_dir,InfoMesh,InfoProblem,InfoMaterial);
    %% find isotherm
    doPlot = 0;
    isotherm_value = InfoProblem.T_LAB;
    addpath '01_plots'
    interphasePoints2.LABx = interphasePoints(1,:);
    interphasePoints2.LABy = interphasePoints(2,:);
    iterations = 0;
    varargin.title = ['Temperature iteration ',num2str(iterations)];
    iso_Temp = plot_results(isotherm_value,doPlot,Temp,13,InfoMesh,interphasePoints2,varargin);
    daspect([1 0.5 1])
    if doPlot == 1
        saveas(gcf,['99_Figures_Afonso\3_Temp1','.pdf'])
    end 
    %% begin loop
    fig_num = 4;
    %nOfPoints_ele = 2;
    addpath '06_OptimizationProcedure'
    plot_points = 0;
    errorT_inGamma = compute_errorOnLAB(Temp,InfoProblem.T_LAB,InfoMesh.list_cut,plot_points,InfoMesh);
    gamma = 0.90; 
    
    Temp_old = T_est;
    LS_temp = 1; 
    
    while errorT_inGamma > 1e-6%* InfoProblem.T_LAB
        iterations = iterations + 1; 
        if LS_temp == 0
            % relaxation among lines (alpha1)
            [newX,newY] = prueba_alpha_LAB_isotherm(gamma,InfoLAB.LABx,InfoLAB.LABy,flip(iso_Temp(1,:)),flip(iso_Temp(2,:)),maxDepth);
            InfoLAB.LABx = newX; 
            InfoLAB.LABy = newY; 
            % compute level sets
            varargin.plotLS = 0; 
            InfoMesh.LS_val_Stokes = LevelSet(InfoMesh.X_v,InfoLAB,varargin);    
            InfoMesh.LS_val_Stokes_unitless = InfoMesh.LS_val_Stokes / InfoProblem.L_ref;
            InfoMesh.LS_val = convertQuad2Lin(InfoMesh.LS_val_Stokes, InfoMesh.X, InfoMesh.T, InfoMesh.T_v);
            daspect([1 1 1])
            saveas(gcf,['99_Figures_Afonso\',num2str(fig_num),'_LS_',num2str(iterations),'.pdf'])
        else
            LS_old = InfoMesh.LS_new; 
            InfoMesh.LS_new = gamma * LS_old + (1-gamma) * (InfoProblem.T_LAB - Temp);
            [InfoMesh.list_Omega1,InfoMesh.list_Omega2,InfoMesh.list_cut, InfoMesh.list_edge_elem1, InfoMesh.list_edge_elem2, PointsInterface] = findElemCrossedLS(tol,plot1,InfoMesh, [] ,InfoProblem);
            daspect([1 1 1])
            PointsInterface(2,:) = maxDepth - PointsInterface(2,:);
            example_plots.LAB = PointsInterface'/InfoProblem.L_ref;
        end 
        % LS- Stokes 
        example_plots.fig_velo=1;
        % relax temperature
        T_relaxed = gamma * Temp_old + (1-gamma) * Temp; 
        % solve Stokes
        [u_mant,p_omega] = ComputeStokesProblem(T_relaxed,p_omega, InfoMesh,InfoProblem,InfoMaterial,example_plots);
        if example_plots.fig_velo == 1
            daspect([1 1 1])
            saveas(gcf,['99_Figures_Afonso\',num2str(fig_num+1),'_u_conv',num2str(iterations),'.pdf'])
        end 
        % compute temp: 1) re-define parameters and cond. ; 2) compute elements crossed by LAB; 3) convert linear u, 4) solve
        grad_T1 = (InfoProblem.T_LAB-InfoProblem.T_sup)/(maxDepth - mean(PointsInterface(2,:) ) ); % [K/m]
        InfoProblem.k2 = InfoProblem.k1 * grad_T1 / InfoProblem.grad_apprT2; % [W/(m K)]
        if cond_dir == 1
            InfoProblem.T_inf = InfoProblem.T_LAB + InfoProblem.grad_apprT2 * (maxDepth-mean(PointsInterface(2,:)) );
        else
            InfoProblem.q2 = InfoProblem.k2 * InfoProblem.grad_apprT2; % [W/(m^2)]
        end 
        % solve
        Temp = solveTemperatureProblem(u_mant,T_relaxed, p_omega,cond_dir,InfoMesh,InfoProblem,InfoMaterial);
        doPlot = 1;
        % plot results
        InterphaseStruct.LABx = PointsInterface(1,:); 
        InterphaseStruct.LABy = PointsInterface(2,:); 
        varargin.title = ['Temperature iteration ',num2str(iterations)];
        iso_Temp = plot_results(isotherm_value,doPlot,Temp,13,InfoMesh,InterphaseStruct,varargin);
        if doPlot == 1
            daspect([1 1 1])
            saveas(gcf,['99_Figures_Afonso\',num2str(fig_num+2),'_Temp',num2str(iterations),'.pdf'])
        end 
        iso_Temp(2,:) = maxDepth - iso_Temp(2,:); % this is due to plot in depth instead of Y
        fig_num = fig_num+3;
        % compute points on LAB to find temperature on it
        errorT_inGamma = compute_errorOnLAB(Temp,InfoProblem.T_LAB,InfoMesh.list_cut,plot_points,InfoMesh);
    end
    
    figure(120); 
    saveas(gcf,['99_Figures_Afonso\nel_x_',num2str(InfoMesh.nel_x),'.pdf'])
    InfoLAB.LABx = flip(iso_Temp(1,:));
    InfoLAB.LABy = flip(iso_Temp(2,:));

end 