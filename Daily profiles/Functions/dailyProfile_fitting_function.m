%% dailyProfile_fitting_function

% This script is developped in relation with the publication 'Modelling E.
% coli removal during real domestic wastewater treatment in outdoor pilot
% scale High Rate Algal Ponds' with the objective to fit the E. coli decay
% model developped in the publication to the data from daily profiles.

% Briefly, this function imports the input variables and values for fitted
% parameters, computes the predicted E. coli cell counts for the data base
% hence forth imported and computes the SSR and MRAE associated to the
% prediction. The outputs of the function are the SSR, the MRAE, and the
% modelled E. coli cell count (log transformed and not).

%% Function

function [SSR,MRAE,coli_model,log_coli_fit] = dailyProfile_fitting_function(fitted_par,...
    time_data,sun_data,pH_data,temp_data,coli_data,n_data,...
    TSS,Q_IN,C_IN,slope_sigma,intercept_sigma,d,S,...
    n_exp)

    k_pH = fitted_par(1);
    teta_pH  = fitted_par(2);
    k_dark = fitted_par(3);
    teta_dark = fitted_par(4);
    alpha_sun = fitted_par(5);

    coli_model = cell(n_exp,1);


    coli_fit = [];
    coli_fit_no1st = [];
    for i_model = 1:n_exp
        coli_model{i_model} = NaN(n_data{i_model},1);

        k = 1; % This loop looks for the first value measured of E. coli. Because the data was formatted, it should not be necessary and it should stop at 1 at all times
        while isnan(coli_data{i_model}(k)) && k < n_data{i_model}
            k = k + 1;
        end
        coli_model{i_model}(k) = coli_data{i_model}(k);

        for m = k + 1:n_data{i_model}
            k_dark_inst = k_dark*teta_dark^(temp_data{i_model}(m-1)-20);
            k_pH_inst = k_pH*teta_pH^(temp_data{i_model}(m-1)-20)*...
                10^(pH_data{i_model}(m-1)-14);
            k_sun_inst = alpha_sun*sun_data{i_model}(m-1)/((slope_sigma*TSS(i_model)+intercept_sigma)*d)*(1-exp(-(slope_sigma*TSS(i_model)+intercept_sigma)*d));

            coli_model{i_model}(m) = coli_model{i_model}(m-1) + ((Q_IN(i_model)/(S*d))*(C_IN(i_model) - coli_model{i_model}(m-1)) - coli_model{i_model}(m-1)*...
                (k_dark_inst + k_pH_inst + k_sun_inst))*(time_data{i_model}(m)-time_data{i_model}(m-1));
            if coli_model{i_model}(m) < 1000
                n_subdiv = 500;
                time_local = time_data{i_model}(m-1):(time_data{i_model}(m) - time_data{i_model}(m-1))/n_subdiv:time_data{i_model}(m);
                temp_local = interp1([time_data{i_model}(m-1),time_data{i_model}(m)],[temp_data{i_model}(m-1),temp_data{i_model}(m)],time_local);
                pH_local = interp1([time_data{i_model}(m-1),time_data{i_model}(m)],[pH_data{i_model}(m-1),pH_data{i_model}(m)],time_local);
                sun_local = interp1([time_data{i_model}(m-1),time_data{i_model}(m)],[sun_data{i_model}(m-1),sun_data{i_model}(m)],time_local);
                coli_model_local = NaN(n_subdiv,1);
                coli_model_local(1) = coli_model{i_model}(m-1);
                for i_subdiv = 2:n_subdiv +1
                    k_dark_inst = k_dark*teta_dark^(temp_local(i_subdiv)-20);
                    k_pH_inst = k_pH*teta_pH^(temp_local(i_subdiv)-20)*10^(pH_local(i_subdiv)-14);
                    k_sun_inst = alpha_sun*sun_local(i_subdiv)/((slope_sigma*TSS(i_model)+intercept_sigma)*d)*(1-exp(-(slope_sigma*TSS(i_model)+intercept_sigma)*d));
                    % Determination of E. coli cell count
                    coli_model_local(i_subdiv) = coli_model_local(i_subdiv-1) + ((Q_IN(i_model)/(d*S))*(C_IN(i_model) - coli_model_local(i_subdiv-1))...
                        - coli_model_local(i_subdiv-1)*(k_dark_inst + k_pH_inst + k_sun_inst))*(time_local(i_subdiv)-time_local(i_subdiv-1));
                end
                coli_model{i_model}(m) = coli_model_local(n_subdiv+1);
            end
        end
        index_coli = find(~isnan(coli_data{i_model}));
        coli_fit = [coli_fit ; coli_data{i_model}(index_coli), coli_model{i_model}(index_coli), i_model*ones(length(index_coli),1)];
        coli_fit_no1st = [coli_fit_no1st;[NaN ; coli_model{i_model}(index_coli(2:end))]]; 
    end
    log_coli_fit = log10(coli_fit);
    log_coli_data = log_coli_fit(:,1);
    log_coli_fit_no1st = log10(coli_fit_no1st);


    SSR = sum((log_coli_fit(:,1)-log_coli_fit(:,2)).^2);
    MRAE = mean(...
        abs(log_coli_data(isnan(log_coli_fit_no1st) == 0) - log_coli_fit_no1st(isnan(log_coli_fit_no1st) == 0))...
        ./log_coli_data(isnan(log_coli_fit_no1st) == 0));

end
    
    