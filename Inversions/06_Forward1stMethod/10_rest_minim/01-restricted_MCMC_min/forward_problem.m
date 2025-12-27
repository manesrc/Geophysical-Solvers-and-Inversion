function [Likelihood,misfit] = forward_problem(g_trial,flux_obs,std_dev,B_0,v_0)

data_model = B_0*g_trial - v_0;
differences = 0;
for ii  = 1:size(flux_obs,2)
    dif_ii = data_model - flux_obs(:,ii);
    differences = differences + dif_ii' * dif_ii;
end 
misfit = (1/(2*std_dev^2)) * differences;
Likelihood = exp(-misfit);
