function [likelihood_vect,misfit_vect,acept_perc,due2prior_reject,due2like_reject] = run_MCMC_minimization(nsimu,nDOFQ,prior_g,g0,n_obs, mean_obs, sigma_obs, B_0,v_0,D_mat,D_matQ, MAT2_new2, const_vect, DOF2, InfoProblem, InfoLAB,InfoMesh, ii_unr, structure)
    % Restricted minimization using MCMC 
    pert_matrix = zeros(nDOFQ,nsimu);   % save each chain obtained
    likelihood_vect = zeros(1,nsimu);           % evolution of likelihood with steps
    misfit_vect = zeros(size(likelihood_vect)); % evolution of the misfit with steps
    %prior_g = [Lower_bound  Upper_bound];   % where the result of g has to be 
    
    %g1 = (Lower_bound + Upper_bound)/2;  % initial guess
    pert_matrix(:,1) = g0; % save in matrix
    sigma = diag(prior_g(:,2)-prior_g(:,1))/10;    % range of possible values for g 
    %sigma = eye(nDOFQ) * (Upper_bound(1)-Lower_bound(1));
    data_obs = mean_obs + sigma_obs*randn(size(B_0,1),n_obs); % synthetic observation data  [ mean 0 variance of 0.05]
    
    dimensionalize_flux = (InfoProblem.T_ref * InfoProblem.k_ref) / InfoProblem.L_ref;
    dimensionalize_gradient2 = dimensionalize_flux /(InfoProblem.k_ref * InfoProblem.k2);

    % calculate initial likelihood
    %[Likelihood,misfit,model_results]
    [likelihood_vect(1),misfit_vect(1)] = forward_problem(g0,data_obs,sigma_obs,B_0,v_0);
    
    
    acept_perc = 0;
    due2prior_reject = 0;
    due2like_reject = 0;
    
    
    for ii = 2:nsimu
        g_ini = pert_matrix(:,ii-1);
        like_ini = likelihood_vect(ii-1);
        new_g = modifyProposedg(sigma,g_ini);
        if any(new_g < prior_g(:,1) | new_g > prior_g(:,2))
            %probPrior = 0;
            pert_matrix(:,ii) = g_ini';
            likelihood_vect(ii) = like_ini;
            misfit_vect(ii) = misfit_vect(ii-1);
            due2prior_reject = due2prior_reject+1;
        else
            [likelihood_g_new,e_m] = forward_problem(new_g,data_obs,sigma_obs,B_0,v_0);
            % which g to save: using log criterion
            log_rm = misfit_vect(ii-1) - e_m;
            if log_rm >= 0
                pert_matrix(:,ii) = new_g';
                likelihood_vect(ii) = likelihood_g_new;
                misfit_vect(ii) = e_m;
                acept_perc = acept_perc + 1;
            else
                pert_matrix(:,ii) = g_ini';
                likelihood_vect(ii) = like_ini;
                misfit_vect(ii) = misfit_vect(ii-1);
                due2like_reject = due2like_reject+1;
            end     
        end 
    end 
    disp(['El % de aceptacion es ',num2str(100*(acept_perc/nsimu))])
    disp(['El % de rechazo por muestrear fuera del rango del prior es ',num2str(100*(due2prior_reject/nsimu))])
    disp(['El % de rechazo por likelihood es ',num2str(100*(due2like_reject/nsimu))])
    
    pert_MCMC = mean(pert_matrix(:,end-(ceil(nsimu/10)):end),2);
    
    
    if structure.Plots ==  1
        line_int_vect = structure.line_int_vect;
        gradient_pertMCMC = (pert_MCMC./line_int_vect) * dimensionalize_gradient2;
        gradient_method2_dim = structure.grad2_minim_analyt_dim;
        x_axis1 = linspace(0,InfoProblem.L_ref/1000,nDOFQ);
        figure(12392); clf;
        plot(x_axis1,gradient_pertMCMC,'r-'); 
        hold on
        plot(x_axis1,gradient_method2_dim,'m--');
        plot(x_axis1,prior_g(:,1)*dimensionalize_gradient2,'k--')
        plot(x_axis1,prior_g(:,2)*dimensionalize_gradient2,'b.-')
        set(gca,'FontSize',16)
        ylabel('pert. to \nabla T [K/km]','FontSize',16)
        xlabel('X [km] ','FontSize',16)
        title('Comparing Matlab rest. min. with MCMC results')
        legend('MCMC','Restricted LS','\delta \nabla T = -0.2 K/km','\delta \nabla T = 0.1 K/km','FontSize',18,'Location','Best')
        axis tight
        grid on
        grid minor
        % plot
        perturbation_longMC = D_mat*pert_MCMC;
        %pertG_MCMC = D_matQ*pert_MCMC;
        vect2_rest_MC = const_vect+ perturbation_longMC;
        % temperatures in Omega2
        T_2_nd_rest_MC = MAT2_new2\vect2_rest_MC;
        % fluxes
        
        matG2 = structure.matG2;
        fluxOmega2_rest_MC = (matG2*T_2_nd_rest_MC) * dimensionalize_flux;
        %plot: 
        ii_r=ii_unr+18;
        results.fluxOmega1_d = structure.flux1;
        results.gradient_T_mean = structure.grad2mean;
        results.T1nd = structure.T_1_nd;
        results.fluxOmega2_gmean_d = structure.flux2mean;
        results.PosMat = structure.PosMat;
        results.DOF1 = structure.DOF1;
        results.DOF2 = DOF2;
        results.gradient_T_mean = structure.grad2mean;
        varargin.restricted = 1;
        plot_temp_flux(ii_r,T_2_nd_rest_MC,pert_MCMC,fluxOmega2_rest_MC,results,InfoMesh,InfoLAB,InfoProblem,varargin)

        if structure.Plots_hist == 1
            n1 = 22;
            mean_DOF_val = zeros(nDOFQ,1);
            var_DOF_val = zeros(nDOFQ,1);
            for kk = 1:nDOFQ
                x_plot = pert_matrix(kk,:);
                mean_DOF_val(kk) = mean(x_plot);
                var_DOF_val(kk) = var(x_plot);
                figure(n1); clf;
                histogram(x_plot)
                n1 = n1+1;
            end
        end 

    end 
