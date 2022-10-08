%% dailyProfile_fitting_function_SSR

% This function enables to compute the SSR directly based on the function
% dailyProfile_fitting_function. This function enables to minimize the SSR
% of the predictive model developped in the publication 'Modelling E.
% coli removal during real domestic wastewater treatment in outdoor pilot
% scale High Rate Algal Ponds' in order to compute the best fit relatively
% to the fitting parameters.

function SSR = dailyProfile_fitting_function_SSR(fitted_par,...
    time_data,sun_data,pH_data,temp_data,coli_data,n_data,...
    TSS,Q_IN,C_IN,slope_sigma,intercept_sigma,d,S,...
    n_exp)

%% This function aims at producing only the SSR output from the dp_fitting_strategy_6 function
% It has the same inputs and basically exectutes this very function.
    [SSR,~] = dailyProfile_fitting_function(fitted_par,...
        time_data,sun_data,pH_data,temp_data,coli_data,n_data,...
        TSS,Q_IN,C_IN,slope_sigma,intercept_sigma,d,S,...
        n_exp);
end


