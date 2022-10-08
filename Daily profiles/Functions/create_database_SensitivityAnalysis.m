%% create_database_SensitivityAnalysis

% This script aims at preparing the inputs for a single sensitivity
% analysis simulation in the objective of fitting uncertainty analysis as
% performed in the files Main_sensitivityAnalysis_fittedParameters and
% Main_sensitivityAnalysis_model_fitness in the context of the publication
% 'Modelling E. coli removal during real domestic wastewater treatment in 
% outdoor pilot scale High Rate Algal Ponds'.

% Briefly, the base value and values of uncertainties are used as input of
% the function. An errored value of each input is calculated based on the
% error used as input.
% NB: the error input is set in the Main scripts calling this function to
% match the needs of the sensitivity analysis.
% The outputs are the transformed values of the input varaibles.

%% Function

function [time_data_MC,sun_data_MC,pH_data_MC,temp_data_MC,coli_data_MC,TSS_MC,slope_sigma_MC,intercept_sigma_MC,d_MC,S_MC,Q_IN_MC,C_IN_MC] =...
    create_database_SensitivityAnalysis(time_data,sun_data,pH_data,temp_data,coli_data,coli_data_IC_minus,coli_data_IC_plus,n_data,...
        TSS_0,C_IN_0,Q_IN_0,slope_sigma_0,intercept_sigma_0,d_0,S_0,...
        err_pH,err_temp,err_sun,err_C_IN,err_coli_HRAP,err_TSS,err_slope_sigma,err_intercept_sigma,err_depth,err_S,err_Q_IN,...
        n_exp)


    time_data_MC = cell(n_exp,1);
    pH_data_MC = cell(n_exp,1);
    temp_data_MC = cell(n_exp,1);
    sun_data_MC = cell(n_exp,1);
    coli_data_MC = cell(n_exp,1);
    
    for i = 1:n_exp
        time_data_MC{i} = time_data{i};
        sun_data_MC{i} = sun_data{i}.*(1+err_sun);
        pH_data_MC{i} = pH_data{i}+err_pH;
        temp_data_MC{i} = temp_data{i} + err_temp;
        if err_coli_HRAP == 1
            coli_data_MC{i} = coli_data_IC_plus{i};
        else
            if err_coli_HRAP == -1
                coli_data_MC{i} = coli_data_IC_minus{i};
            else
                coli_data_MC{i} = coli_data{i} + err_coli_HRAP*(coli_data_IC_plus{i} - coli_data_IC_minus{i})/2;
            end
        end    
        coli_data_MC{i}(coli_data_MC{i} < 0) = 1;
    end
    
    
    TSS_MC = TSS_0.*(1 + err_TSS);
    slope_sigma_MC = slope_sigma_0 + err_slope_sigma;
    intercept_sigma_MC = intercept_sigma_0 + err_intercept_sigma;
    d_MC = d_0 + err_depth;
    S_MC = S_0 + err_S;
    Q_IN_MC = Q_IN_0*(1 + err_Q_IN);
    C_IN_MC = max(1,C_IN_0*(1 + err_C_IN));
    
end
    
    