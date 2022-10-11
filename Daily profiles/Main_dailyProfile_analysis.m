%% Main_dailyProfile_analysis

% This script aims at performing most analysis realeted withdaily profile
% experiments in the publication 'Modelling E. coli removal during real
% domestic wastewater treatment in outdoor pilot scale High Rate Algal
% Ponds'

% This script is hence used to draw Figures 2, S5.1, S5.2, S7.1, S6.2.  

%% Options
clear all
fs = 14;
formatSpec = '%.2f'; % option which limits to 2 digits the results presented in figure titles

%% Figure S5.1 and S5.2
% These figures represent the fit of E. coli cell count during daily
% profile analysis witht he original fit from the publication Chambonniere
% et al. 2022.

    % Import of the data needed

    % Initialization
        % Loading the data to be able to compute E. coli cell count during
        % pilot scale daily profile.

        % Measured parameters

        measuredParametersTable = xlsread('./Daily profiles/Daily_profile_transformedData.xlsx','Model');
        % depth (m)
        d = measuredParametersTable(7,2);
        % Pond surface (m2)
        S = measuredParametersTable(8,2);
        % Parameters from the linear fit between TSS (mg/L) and light
        % attenuation coefficient (m-1)
        slope_sigma = 0.159;
        intercept_sigma = 11.769;
        % TSS concentration (mg/L)
        TSSTable  = measuredParametersTable(10:16,2);
        % Influent flow rate (m3/d)
        QIN = measuredParametersTable(18:24,3);
        % Influent E. coli cell count (MPN/100 mL)
        CIN = measuredParametersTable(26:32,2);
        % Number of daily profile experiments carried out 
        nExp = 7;
        % nomber of E. coli cell count measurement total
        nDataColi = 80; 

        % Measured data import
        measuredVariables = {xlsread('./Daily profiles/Daily_profile_transformedData.xlsx','30_09'),...
            xlsread('./Daily profiles/Daily_profile_transformedData.xlsx','12_10'),...
            xlsread('./Daily profiles/Daily_profile_transformedData.xlsx','29_10'),...
            xlsread('./Daily profiles/Daily_profile_transformedData.xlsx','17_11'),...
            xlsread('./Daily profiles/Daily_profile_transformedData.xlsx','03_02'),...
            xlsread('./Daily profiles/Daily_profile_transformedData.xlsx','10_02'),...
            xlsread('./Daily profiles/Daily_profile_transformedData.xlsx','16_03')};



        % Creation of table for storage of measured  variables
        % nomber of measurement for environmental parameters per experiment
        nData = cell(nExp,1); 
        % Table of time data (d)
        time_data = cell(nExp,1);
        % Table of sunlight direct intensity (W/m2)
        sun_data = cell(nExp,1);
        % Table of dissolved oxygen (mg/L)
        DO_data = cell(nExp,1); 
        % Table of pH
        pH_data = cell(nExp,1);
        % Table of temperature (°C)
        temp_data = cell(nExp,1);
        % Table of E. coli cell count measured (MPN/100mL)
        coli_data = cell(nExp,1);
         % Table of E. coli cell count predicted (MPN/100mL)
        coli_model = cell(nExp,1);
         % Table of E. coli cell count predicted lower CI95% value (MPN/100mL)
        coli_data_IC_minus = cell(nExp,1);
         % Table of E. coli cell count predicted upper CI95% value (MPN/100mL)
        coli_data_IC_plus = cell(nExp,1);

        % Storage of the input data in cell array. Each cell correspond to a
        % daily profile ordered in chronological order. Because E. coli cell
        % count was measured only at certain times from samples, NaN values in
        % the data vector must be filtered out to obtain only the data. These
        % values are also stored in vector as log-transformed values.
        coli_fit = [];
        coli_fit_IC_minus = [];
        coli_fit_IC_plus = [];
        for i = 1:nExp
            nData{i} = size(measuredVariables{i},1);
            time_data{i} = measuredVariables{i}(:,1);
            sun_data{i} = measuredVariables{i}(:,5);
            pH_data{i} = measuredVariables{i}(:,2);
            temp_data{i} = measuredVariables{i}(:,3);
            coli_data{i} = measuredVariables{i}(:,6);
            coli_data_IC_minus{i} = measuredVariables{i}(:,7);
            coli_data_IC_plus{i} = measuredVariables{i}(:,8);
            coli_model{i} = NaN(nData{i},1);

            index_coli = find(~isnan(coli_data{i}));
            coli_fit = [coli_fit;coli_data{i}(index_coli)];
            coli_fit_IC_minus = [coli_fit_IC_minus;coli_data_IC_minus{i}(index_coli)];
            coli_fit_IC_plus = [coli_fit_IC_plus;coli_data_IC_plus{i}(index_coli)];
        end 
        coli_logcounts = log10(coli_fit);
        coli_logcounts_IC_minus = log10(coli_fit_IC_minus);
        coli_logcounts_IC_plus = log10(coli_fit_IC_plus);


    % Bench scale model validation from Chambonniere et al. 2022 original
    % fit

        % Fitting parameters values as obtained from Chambonniere et al. (2022)
        k_pH = 2860;
        teta_pH = 1.45;
        alpha_sun = 0;
        k_dark = 47.4;
        teta_dark = 1.00;

        % In the following, E. coli cells predicted from Chambonniere et al.
        % 2022 fitted parameters and input parameters/variables are computed
        coli_fit = [];
        coli_fit_no1st = []; % This vector will store the predicted E. coli cell count when removing the first value of each day in order to compute the fitness parameters
        for i_model = 1:nExp

            % The first existing data for E. coli cell count during an
            % experiment is looked for to initialize the model calculation
            k = 1;
            while isnan(coli_data{i_model}(k)) && k < nData{i_model}
                k = k + 1;
            end
            coli_model{i_model}(k) = coli_data{i_model}(k);

            % Predicted E. coli cell count are now determined usin Euler method 
            for m = k + 1:nData{i_model}
                % The decay rate induced by each mechanisms are first computed
                % (d-1)
                k_dark_inst = k_dark*teta_dark^(temp_data{i_model}(m-1)-20);
                k_pH_inst = k_pH*teta_pH^(temp_data{i_model}(m-1)-20)*10^(pH_data{i_model}(m-1)-14);
                k_sun_inst = alpha_sun*sun_data{i_model}(m-1)/((slope_sigma*TSSTable(i_model)+intercept_sigma)*d)*(1-exp(-(slope_sigma*TSSTable(i_model)+intercept_sigma)*d));
                % Prediction of E. coli cell count
                coli_model{i_model}(m) = coli_model{i_model}(m-1)*(1 - (k_dark_inst + k_pH_inst + k_sun_inst)*(time_data{i_model}(m)-time_data{i_model}(m-1))) + ...
                    (CIN(i_model)- coli_model{i_model}(m-1))*QIN(i_model)/(S*d)*(time_data{i_model}(m)-time_data{i_model}(m-1));
                % Because the Euler method would sometimes lead to negative
                % predicted E. coli cell counts, it was added this loop in
                % which the step of calculation is reduced when E. coli cell
                % count would get below a certain value (here 1000 MPN/100 mL).
                % The calculation step is divided and all measured data are
                % linearly interpolated on the same step of time.
                % The calculation is then resumed as previously.
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
                        k_sun_inst = alpha_sun*sun_local(i_subdiv)/((slope_sigma*TSSTable(i_model)+intercept_sigma)*d)*(1-exp(-(slope_sigma*TSSTable(i_model)+intercept_sigma)*d));
                        % Determination of E. coli cell count
                        coli_model_local(i_subdiv) = coli_model_local(i_subdiv-1)*(1 - (k_dark_inst + k_pH_inst + k_sun_inst)*(time_local(i_subdiv)-time_local(i_subdiv-1))) + ...
                            (CIN(i_model)- coli_model_local(i_subdiv-1))*QIN(i_model)/(S*d)*(time_local(i_subdiv)-time_local(i_subdiv-1));
                    end
                    coli_model{i_model}(m) = coli_model_local(n_subdiv+1);
                end
            end
            % Storage of the results    
            index_coli = find(~isnan(coli_data{i_model}));
            coli_fit = [coli_fit ; coli_model{i_model}(index_coli)];
            coli_fit_no1st = [coli_fit_no1st;[NaN ; coli_model{i_model}(index_coli(2:end))]]; 

        end

        % Log transformation of the predicted E. coli cell count and
        % calculation of fitness performance parameters (MRAE, R2)
        coli_predicted_logcounts = log10(coli_fit_no1st);
        MRAE = 100*mean(...
            abs(coli_logcounts(isnan(coli_predicted_logcounts) == 0) - coli_predicted_logcounts(isnan(coli_predicted_logcounts) == 0))...
            ./coli_logcounts(isnan(coli_predicted_logcounts) == 0));
        R_squared = 1 - nansum((coli_predicted_logcounts-coli_logcounts).^2)/...
            sum((coli_logcounts(isnan(coli_predicted_logcounts) == 0) - ...
            mean(coli_logcounts(isnan(coli_predicted_logcounts) == 0))).^2);


    % Figures of Chambonniere et al. 2022 model fitness

    % Figure S5.1
        figure(1), clf, hold on
        plot(coli_predicted_logcounts,coli_logcounts,'ob','MarkerSize',6,'LineWidth',2)
        plot([2,3,4,5,6,7],[2,3,4,5,6,7],'-k','LineWidth',2)
        errorbar(coli_predicted_logcounts,coli_logcounts,...
                coli_logcounts_IC_minus - coli_logcounts,...
                coli_logcounts_IC_plus - coli_logcounts,...
                '.k','MarkerSize',0.1,'LineWidth',1)
        ylabel({'Measured {\itE. coli} cell counts' ; '(log_1_0 MPN{\cdot}100 mL^-^1)'},'FontSize',fs,...
            'FontWeight','bold','FontName','Arial')
        xlabel({'Predicted {\itE. coli} cell counts' ; '(log_1_0 MPN{\cdot}100 mL^-^1)'},'FontSize',fs,...
            'FontWeight','bold','FontName','Arial')
        ax2 = gca;
        ax2.FontSize = fs-2;
        ax2.XLim = [2 7];
        ax2.YLim = [2 7];
        title(strcat('R^2 = ',{' '},num2str(R_squared,formatSpec),{' '},';',...
            {' '},'MRAE =',{' '},num2str(MRAE,formatSpec),{' '},'%'),...
            'FontSize',fs,'FontWeight','bold','FontName','Arial')

        fig = gcf;
        fig.Position = [50 50 500 450];

    
    % Study of Monte Carlo from original fit: the fitted parameters
    % predicted by Chambonniere et al. is analyzed in light of the Monte
    % Carlo analysis they performed to determine fitting uncertainty. 
    
        % The results from the authors Monte Carlo analysis are first imported
        load('./Daily profiles/bench_data_fitting_MC_analysis.mat',...
            'k_pH_MC','teta_pH_MC','alpha_sun_MC','k_nat_MC','teta_nat_MC','N_MC')

        % Creation of tables used for storage of the inputs/ouputs for each
        % Monte Carlo loop
        coli_predicted_logcounts = NaN(nDataColi,N_MC);
        coli_model_MC = cell(nExp,1); % Table of E. coli cell count predicted (MPN/100mL)
        for i = 1:nExp
            coli_model_MC{i} = NaN(nData{i},N_MC);
        end

        % Loop for the calculation of predicted E. coli cell count with each
        % set of parameters from Chambonniere et al. 2022 Monte Carlo run
        for i_MC = 1:N_MC

            k_pH = k_pH_MC(i_MC);
            teta_pH = teta_pH_MC(i_MC);
            alpha_sun = alpha_sun_MC(i_MC);
            k_dark = k_nat_MC(i_MC);
            teta_dark = teta_nat_MC(i_MC);

            coli_fit = [];
            coli_fit_no1st = [];
            for i_model = 1:nExp

                k = 1;
                while isnan(coli_data{i_model}(k)) && k < nData{i_model}
                    k = k + 1;
                end
                coli_model_MC{i_model}(k,i_MC) = coli_data{i_model}(k);


                for m = k + 1:nData{i_model}
                    % Determination of per mechanism decay rate
                    k_dark_inst = k_dark*teta_dark^(temp_data{i_model}(m-1)-20);
                    k_pH_inst = k_pH*teta_pH^(temp_data{i_model}(m-1)-20)*10^(pH_data{i_model}(m-1)-14);
                    k_sun_inst = alpha_sun*sun_data{i_model}(m-1)/((slope_sigma*TSSTable(i_model)+intercept_sigma)*d)*(1-exp(-(slope_sigma*TSSTable(i_model)+intercept_sigma)*d));
                    % Determination of E. coli cell count
                    coli_model_MC{i_model}(m,i_MC) = coli_model_MC{i_model}(m-1,i_MC)*(1 - (k_dark_inst + k_pH_inst + k_sun_inst)*(time_data{i_model}(m)-time_data{i_model}(m-1))) + ...
                        (CIN(i_model)- coli_model_MC{i_model}(m-1,i_MC))*QIN(i_model)/(S*d)*(time_data{i_model}(m)-time_data{i_model}(m-1));
                    n_subdiv = 250;
                    if coli_model_MC{i_model}(m,i_MC) < 1000
                        time_local = time_data{i_model}(m-1):(time_data{i_model}(m) - time_data{i_model}(m-1))/n_subdiv:time_data{i_model}(m);
                        temp_local = interp1([time_data{i_model}(m-1),time_data{i_model}(m)],[temp_data{i_model}(m-1),temp_data{i_model}(m)],time_local);
                        pH_local = interp1([time_data{i_model}(m-1),time_data{i_model}(m)],[pH_data{i_model}(m-1),pH_data{i_model}(m)],time_local);
                        sun_local = interp1([time_data{i_model}(m-1),time_data{i_model}(m)],[sun_data{i_model}(m-1),sun_data{i_model}(m)],time_local);
                        coli_model_local = NaN(n_subdiv,1);
                        coli_model_local(1) = coli_model_MC{i_model}(m-1,i_MC);
                        for i_subdiv = 2:n_subdiv +1
                            k_dark_inst = k_dark*teta_dark^(temp_local(i_subdiv)-20);
                            k_pH_inst = k_pH*teta_pH^(temp_local(i_subdiv)-20)*10^(pH_local(i_subdiv)-14);
                            k_sun_inst = alpha_sun*sun_local(i_subdiv)/((slope_sigma*TSSTable(i_model)+intercept_sigma)*d)*(1-exp(-(slope_sigma*TSSTable(i_model)+intercept_sigma)*d));
                            % Determination of E. coli cell count
                            coli_model_local(i_subdiv) = coli_model_local(i_subdiv-1)*(1 - (k_dark_inst + k_pH_inst + k_sun_inst)*(time_local(i_subdiv)-time_local(i_subdiv-1))) + ...
                                (CIN(i_model)- coli_model_local(i_subdiv-1))*QIN(i_model)/(S*d)*(time_local(i_subdiv)-time_local(i_subdiv-1));
                        end
                    coli_model_MC{i_model}(m,i_MC) = coli_model_local(n_subdiv+1);
                    end
                end

                index_coli = find(~isnan(coli_data{i_model}));
                coli_fit = [coli_fit ; coli_model_MC{i_model}(index_coli,i_MC)];
                coli_fit_no1st = [coli_fit_no1st;[NaN ; coli_model_MC{i_model}(index_coli(2:end),i_MC)]]; 

            end

            coli_predicted_logcounts(:,i_MC) = log10(coli_fit_no1st);

        end

        % Uncertainty fit: the following data table are created to show certain
        % percentiles of the values predicted from all the predicted values
        % with Monte Carlo generated sets of fitting parameters.

        median_curve = cell(nExp,1);
        area_table = cell(nExp,1);
        prctile25_curve = cell(nExp,1);
        prctile75_curve = cell(nExp,1);
        prctile5_curve = cell(nExp,1);
        prctile95_curve = cell(nExp,1);
        max_curve = cell(nExp,1);
        min_curve = cell(nExp,1);
        for iExp = 1:nExp
            median_curve{iExp} = NaN(size(coli_model_MC{iExp},1),1);
            prctile25_curve{iExp} = NaN(size(coli_model_MC{iExp},1),1); 
            prctile75_curve{iExp} = NaN(size(coli_model_MC{iExp},1),1); 
            prctile5_curve{iExp} = NaN(size(coli_model_MC{iExp},1),1); 
            prctile95_curve{iExp} = NaN(size(coli_model_MC{iExp},1),1); 
            max_curve{iExp} = NaN(size(coli_model_MC{iExp},1),1);
            min_curve{iExp} = NaN(size(coli_model_MC{iExp},1),1);
            for iPoint = 1:size(coli_model_MC{iExp},1)
                median_curve{iExp}(iPoint) = prctile(log10(coli_model_MC{iExp}(iPoint,:)),50);
                prctile25_curve{iExp}(iPoint) = prctile(log10(coli_model_MC{iExp}(iPoint,:)),25);
                prctile75_curve{iExp}(iPoint) = prctile(log10(coli_model_MC{iExp}(iPoint,:)),75);
                prctile5_curve{iExp}(iPoint) = prctile(log10(coli_model_MC{iExp}(iPoint,:)),5);
                prctile95_curve{iExp}(iPoint) = prctile(log10(coli_model_MC{iExp}(iPoint,:)),95);
                max_curve{iExp}(iPoint) = max(log10(coli_model_MC{iExp}(iPoint,:)));
                min_curve{iExp}(iPoint) = min(log10(coli_model_MC{iExp}(iPoint,:)));
            end
            area_table{iExp} = [ min_curve{iExp} ,...
                prctile5_curve{iExp} - min_curve{iExp},...
                prctile25_curve{iExp} - prctile5_curve{iExp}, ...
                median_curve{iExp} - prctile25_curve{iExp},...
                prctile75_curve{iExp} - median_curve{iExp} ,...
                prctile95_curve{iExp} - prctile75_curve{iExp},...
                max_curve{iExp} - prctile95_curve{iExp}];
        end

        % Figure S5.2 
        figure(3), clf, hold on
        colours_MC = [1 1 1].*[1 ; 0.95 ; 0.75 ; 0.5 ; 0.5 ; 0.75; 0.95];
        for i = 1:nExp
            subplot(2,4,i), hold on
            time_data_loc = datetime(0,1,floor(time_data{i}),0,0,3600*24*(time_data{i} -floor(time_data{i})));
            aplot = area(datenum(time_data_loc),area_table{i},1,'ShowBaseLine','off');
            for iColor = 1:7
                aplot(iColor).FaceColor = colours_MC(iColor,:);
                aplot(iColor).LineStyle = 'none';
            end
            bplot = plot(datenum(time_data_loc),median_curve{i},'-k','LineWidth',2) ;
            cplot = plot(datenum(time_data_loc),log10(coli_data{i}),'ob','MarkerSize',6,'LineWidth',2) ;
            dplot = errorbar(datenum(time_data_loc),log10(coli_data{i}),...
                log10(coli_data{i}) - log10(coli_data_IC_minus{i}),...
                log10(coli_data_IC_plus{i}) - log10(coli_data{i}),...
                '.k','MarkerSize',0.1,'LineWidth',1);
            ylim([0 7])
            datetick('x')
            ax1 = gca;
            ax1.FontSize = fs-4;
            ax1.XLim = [datenum(dateshift(time_data_loc(1),'start','day') + 8/24),...
                datenum(dateshift(time_data_loc(1),'start','day') + 1 + 10/24)];
            xticks(ax1,[datenum(dateshift(time_data_loc(1),'start','day') + 9/24),...
                datenum(dateshift(time_data_loc(1),'start','day') + 21/24), ...
                datenum(dateshift(time_data_loc(1),'start','day') + 1 + 9/24)]);
            xticklabels(ax1,{datestr(dateshift(time_data_loc(1),'start','day') + 9/24,'dd/mm HH:MM'),...
                datestr(dateshift(time_data_loc(1),'start','day') + 21/24,'dd/mm HH:MM'), ...
                datestr(dateshift(time_data_loc(1),'start','day') + 1 + 9/24,'dd/mm HH:MM')});
            xtickangle(ax1,45)
            if (i >= 1) && (i <= 4)
                title('Spring','FontSize',fs - 2,'FontWeight','bold')
            else
                if (i >= 5) && (i <= 6)
                    title('Summer','FontSize',fs - 2,'FontWeight','bold')
                else
                    title('Fall','FontSize',fs - 2,'FontWeight','bold')
                end
            end
            set(ax1,'Box','off')
        end
        legend1 = legend([aplot(2) aplot(3) aplot(4) bplot(1) cplot(1)],{'Outside 5/95 percentile','5-25 or 75-95 percentile',...
            '25 - 75 percentile','Median','Measured data'},'FontSize',fs-2,'FontWeight','bold',...
            'FontName','Arial')
        set(legend1,...
            'Position',[0.740103715986764 0.238425922487196 0.168619794895252 0.152243593182319]);
        subplot(2,4,1)
        ylabel({'{\itE. coli} cell counts' ; '(log_1_0 MPN{\cdot}100 mL^-^1)'},'FontSize',fs,...
            'FontWeight','bold','FontName','Arial')
        subplot(2,4,5)
        ylabel({'{\itE. coli} cell counts' ; '(log_1_0 MPN{\cdot}100 mL^-^1)'},'FontSize',fs,...
            'FontWeight','bold','FontName','Arial')

        fig = gcf;
        fig.Position = [50 50 1100 700];


%% Part II: Fit of data based on only 3 daily profiles: 
% In this part, the new calibration of the model is performed based on the
% 3 daily profile of the fitting data set.

% Name of the fitting data sets excel sheet to limit the data import
dp_name = {'30_09' , '29_10' , '03_02'};
% TSS concentration values (mg/L)
TSSTable = [248 , 295 , 265 ];
% Inlet E. coli by experiment (MPN/100 mL)
CIN = [4.36 , 3.94 , 26 ]*10^6;
% Inlet flowrate by experiment (m3/d)
QIN = [0.1056 , 0.0984 , 0.144 ];
% depth (m)
d = 0.25;
% Surface (m2)
S = 3.42;
% Number of experiment used for this analysis
nExp = 3;
slope_sigma = 0.159;
intercept_sigma = 11.769;
        
% Storage of input data: the input data vectors are stored in cell array,
% each dimension corresponding to a given day.
time_data = cell(nExp,1);
pH_data = cell(nExp,1);
temp_data = cell(nExp,1);
sun_data = cell(nExp,1);
coli_data = cell(nExp,1);
nData = cell(nExp,1);
coli_data_IC_minus = cell(nExp,1);
coli_data_IC_plus = cell(nExp,1);

for i = 1:nExp
    A = xlsread('.//Daily profiles/Daily_profile_transformedData.xlsx',dp_name{i});
    
    time_data{i} = A(:,1);
    sun_data{i} = A(:,5);
    pH_data{i} = A(:,2);
    temp_data{i} = A(:,3);
    coli_data{i} = A(:,6);
    coli_data_IC_minus{i} = A(:,7);
    coli_data_IC_plus{i} = A(:,8);
    nData{i} = length(time_data{i});
end

% Fitting algorithm: the fitting algorithm aims at maximizing the sum of
% squared residuals (SSR) when varying the fitting parameters

% The fitted parameters are initialized based on the values calculated in
% the best fit at bench scale of Chambonniere et al. 2022 publication

k_pH_0 = 2860;
teta_pH_0 = 1.45;
k_dark_0 = 47.4;
teta_dark_0 = 1.00;
alpha_sun_0 = 0;

x0 = [k_pH_0,teta_pH_0,k_dark_0,teta_dark_0,alpha_sun_0];

% Parameters for fitting algorithm: the fitting algorithm is performed
% assuming realistic lower [value stored in lb] and upper boundaries [value
% stored in ub] for each fitting parameters.

lb = [0,1,0,1,0];
ub = [500000,10,200,10,10];

% The maximum number of iteration of the algorithm is set at 200, other
% fitting routine parameters are se by default (see Matlab documentation).
A = []; b = []; Aeq = []; beq = []; nonlcon = [];
options = optimoptions('fmincon','MaxIterations',200,'PlotFcn','optimplotfval');

% The fitting routine is performed minimizing the function
% dailyProfile_fitting_function_SSR developped for this study and computing
% the model SSR based on data from the experiments carried out for this
% study.
[xfit,SSR,exitflag,output] = fmincon(@(x)dailyProfile_fitting_function_SSR(x,...
    time_data,sun_data,pH_data,temp_data,coli_data,nData,...
    TSSTable,QIN,CIN,slope_sigma,intercept_sigma,d,S,...
    nExp),...
    x0,...
    A,b,Aeq,beq,lb,ub,nonlcon,...
    options);

% The model fitness optimized is stored.
SSR_tot = SSR;

% The newly determined optimized fitting parameters are stored
k_pH = xfit(1);
teta_pH = xfit(2);
k_dark = xfit(3);
teta_dark = xfit(4);
alpha_sun = xfit(5);

k_pH_bestFit0 = xfit(1);
teta_pH_bestFit0 = xfit(2);
k_dark_bestFit0 = xfit(3);
teta_dark_bestFit0 = xfit(4);
alpha_sun_bestFit0 = xfit(5);


% Based on the fitted data, E. coli cell count as predicted with the
% optimized set of fitting parameters are recomputed in order to draw the
% figures needed. The same method as used prior when assessing CHambonniere
% et al. (2022) model fitness is used.

coli_fit = [];
coli_fit_IC_minus = [];
coli_fit_IC_plus = [];
for i = 1:nExp
    coli_model{i} = NaN(nData{i},1);    
    index_coli = find(~isnan(coli_data{i}));
    coli_fit = [coli_fit;coli_data{i}(index_coli)];
    coli_fit_IC_minus = [coli_fit_IC_minus;coli_data_IC_minus{i}(index_coli)];
    coli_fit_IC_plus = [coli_fit_IC_plus;coli_data_IC_plus{i}(index_coli)];
end 
coli_logcounts = log10(coli_fit);
coli_logcounts_IC_minus = log10(coli_fit_IC_minus);
coli_logcounts_IC_plus = log10(coli_fit_IC_plus);
coli_fit_no1st = [];

for i_model = 1:nExp

    % Finding first existing data of an experiment to initialize the
    % model calculation
    k = 1;
    while isnan(coli_data{i_model}(k)) && k < nData{i_model}
        k = k + 1;
    end
    coli_model{i_model}(k) = coli_data{i_model}(k);

    for m = k + 1:nData{i_model}
        % Determination of per mechanism decay rate
        k_dark_inst = k_dark*teta_dark^(temp_data{i_model}(m-1)-20);
        k_pH_inst = k_pH*teta_pH^(temp_data{i_model}(m-1)-20)*...
            10^(pH_data{i_model}(m-1)-14);
        k_sun_inst = alpha_sun*sun_data{i_model}(m-1)/...
            ((slope_sigma*TSSTable(i_model) + intercept_sigma)*d)*...
            (1-exp(-(slope_sigma*TSSTable(i_model) + intercept_sigma)*d));
        % Determination of E. coli cell count
        coli_model{i_model}(m) = coli_model{i_model}(m-1)*...
            (1 - (k_dark_inst + k_pH_inst + k_sun_inst)*...
            (time_data{i_model}(m)-time_data{i_model}(m-1))) + ...
            (CIN(i_model)- coli_model{i_model}(m-1))*QIN(i_model)/(S*d)*...
            (time_data{i_model}(m)-time_data{i_model}(m-1));
        % Safety loop in case the calculation step is too large
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
                k_pH_inst = k_pH*teta_pH^(temp_local(i_subdiv)-20)*...
                    10^(pH_local(i_subdiv)-14);
                k_sun_inst = alpha_sun*sun_local(i_subdiv)/...
                    ((slope_sigma*TSSTable(i_model) + intercept_sigma)*d)*...
                    (1-exp(-(slope_sigma*TSSTable(i_model) + intercept_sigma)*d));
                % Determination of E. coli cell count
                coli_model_local(i_subdiv) = coli_model_local(i_subdiv-1)*...
                    (1 - (k_dark_inst + k_pH_inst + k_sun_inst)*...
                    (time_local(i_subdiv)-time_local(i_subdiv-1))) + ...
                    (CIN(i_model)- coli_model_local(i_subdiv-1))*...
                    QIN(i_model)/(S*d)*(time_local(i_subdiv)-time_local(i_subdiv-1));
            end
        coli_model{i_model}(m) = coli_model_local(n_subdiv+1);
        end
    end

    index_coli = find(~isnan(coli_data{i_model}));
    coli_fit = [coli_fit ; coli_model{i_model}(index_coli)];
    coli_fit_no1st = [coli_fit_no1st;[NaN ; coli_model{i_model}(index_coli(2:end))]]; 

end

% The fitted parameters are stored in a Matlab data set

save(strcat('./Daily profiles/','fittingParameters_bestFit0.mat'),'k_pH_bestFit0','teta_pH_bestFit0','k_dark_bestFit0','teta_dark_bestFit0','alpha_sun_bestFit0','SSR')

% Figure S6.1 is now drawn
figure(20), clf, hold on
for i = 1:nExp
    subplot(2,2,i), hold on
    time_data_loc = datetime(0,1,floor(time_data{i}),0,0,3600*24*(time_data{i} -floor(time_data{i})));
    plot(datenum(time_data_loc),log10(coli_model{i}),'-k','LineWidth',2)
    plot(datenum(time_data_loc),log10(coli_data{i}),'ob','MarkerSize',6,'LineWidth',2)
    errorbar(datenum(time_data_loc),log10(coli_data{i}),...
        log10(coli_data{i}) - log10(coli_data_IC_minus{i}),...
        log10(coli_data_IC_plus{i}) - log10(coli_data{i}),...
        '.k','MarkerSize',0.1,'LineWidth',1)
    ylim([0 6])
    datetick('x')
    ax1 = gca
    ax1.FontSize = fs-4
    ax1.XLim = [datenum(dateshift(time_data_loc(1),'start','day') + 8/24),...
        datenum(dateshift(time_data_loc(1),'start','day') + 1 + 10/24)];
    xticks(ax1,[datenum(dateshift(time_data_loc(1),'start','day') + 9/24),...
        datenum(dateshift(time_data_loc(1),'start','day') + 21/24), ...
        datenum(dateshift(time_data_loc(1),'start','day') + 1 + 9/24)]);
    xticklabels(ax1,{datestr(dateshift(time_data_loc(1),'start','day') + 9/24,'dd/mm HH:MM'),...
        datestr(dateshift(time_data_loc(1),'start','day') + 21/24,'dd/mm HH:MM'), ...
        datestr(dateshift(time_data_loc(1),'start','day') + 1 + 9/24,'dd/mm HH:MM')});
    xtickangle(ax1,45)
end
subplot(2,2,1)
ylabel({'{\itE. coli} cell counts' ; '(log_1_0 MPN{\cdot}100 mL^-^1)'},'FontSize',fs,...
    'FontWeight','bold','FontName','Arial')
subplot(2,2,3)
ylabel({'{\itE. coli} cell counts' ; '(log_1_0 MPN{\cdot}100 mL^-^1)'},'FontSize',fs,...
    'FontWeight','bold','FontName','Arial')


coli_predicted_logcounts = log10(coli_fit_no1st);
MRAE = 100*mean(...
    abs(coli_logcounts(isnan(coli_predicted_logcounts) == 0) - coli_predicted_logcounts(isnan(coli_predicted_logcounts) == 0))...
    ./coli_logcounts(isnan(coli_predicted_logcounts) == 0));
R_squared = 1 - nansum((coli_predicted_logcounts-coli_logcounts).^2)/...
    sum((coli_logcounts(isnan(coli_predicted_logcounts) == 0) - ...
    mean(coli_logcounts(isnan(coli_predicted_logcounts) == 0))).^2);
% figure(10), clf, hold on
subplot(2,2,4), hold on
plot(coli_predicted_logcounts,coli_logcounts,'ob','MarkerSize',6,'LineWidth',2)
xlim([2.5 6.5]), ylim([2.5 6.5])
plot([2,3,4,5,6,7],[2,3,4,5,6,7],'-k','LineWidth',2)
errorbar(coli_predicted_logcounts,coli_logcounts,...
        coli_logcounts_IC_minus - coli_logcounts,...
        coli_logcounts_IC_plus - coli_logcounts,...
        '.k','MarkerSize',0.1,'LineWidth',1)
ylabel({'Measured {\itE. coli} cell counts' ; '(log_1_0 MPN{\cdot}100 mL)'},'FontSize',fs-2,...
    'FontWeight','bold','FontName','Arial')
xlabel({'Predicted {\itE. coli} cell counts' ; '(log_1_0 MPN{\cdot}100 mL^-^1)'},'FontSize',fs-2,...
    'FontWeight','bold','FontName','Arial')
ax2 = gca;
ax2.FontSize = fs-2;
ax2.XLim = [2 7];
ax2.YLim = [2 7];
title(strcat('R^2 = ',{' '},num2str(R_squared,formatSpec),{' '},';',...
    {' '},'MRAE =',{' '},num2str(MRAE,formatSpec),{' '},'%'),...
    'FontSize',fs-2,'FontWeight','bold','FontName','Arial')

fig = gcf;
fig.Position = [50 50 750 700];


%% Part III: Validation of the fitted model on the remaining 4 daily profiles
% The new set of fitted parameters is now tested to predict E. coli cell
% count on the validation data set. The calculation is performed as during
% prior analysis. This section necessitates previous section to have been
% ran to initiate the values of fitted paramereters from the best fit
% values.

fs = 14;
formatSpec = '%.2f';

dp_name = {'12_10' , '17_11' , '10_02' , '16_03' };
TSSTable = [220 , 310 , 210 , 740 ];
CIN = [3.47 , 6.15 , 4.96 , 3.45 ]*10^6;
QIN = [0.0984 , 0.1056 , 0.1488 , 0.12 ];
nExp = 4;

% Creation of table of measured  variables
time_data = cell(nExp,1);  
sun_data = cell(nExp,1); 
DO_data = cell(nExp,1); 
pH_data = cell(nExp,1); 
temp_data = cell(nExp,1);
coli_data = cell(nExp,1);
coli_model = cell(nExp,1);
coli_data_IC_minus = cell(nExp,1);
coli_data_IC_plus = cell(nExp,1);

coli_fit = [];
coli_fit_IC_minus = [];
coli_fit_IC_plus = [];
for i = 1:nExp
    A = xlsread('./Daily profiles/Daily_profile_transformedData.xlsx',dp_name{i});
    
    time_data{i} = A(:,1);
    sun_data{i} = A(:,5);
    pH_data{i} = A(:,2);
    temp_data{i} = A(:,3);
    coli_data{i} = A(:,6);
    coli_data_IC_minus{i} = A(:,7);
    coli_data_IC_plus{i} = A(:,8);
    nData{i} = length(time_data{i});
    coli_model{i} = NaN(nData{i},1);
    
    index_coli = find(~isnan(coli_data{i}));
    coli_fit = [coli_fit;coli_data{i}(index_coli)];
    coli_fit_IC_minus = [coli_fit_IC_minus;coli_data_IC_minus{i}(index_coli)];
    coli_fit_IC_plus = [coli_fit_IC_plus;coli_data_IC_plus{i}(index_coli)];
end 
coli_logcounts = log10(coli_fit);
coli_logcounts_IC_minus = log10(coli_fit_IC_minus);
coli_logcounts_IC_plus = log10(coli_fit_IC_plus);

coli_fit_no1st = [];
for i_model = 1:nExp

    % Finding first existing data of an experiment to initialize the
    % model calculation
    k = 1;
    while isnan(coli_data{i_model}(k)) && k < nData{i_model}
        k = k + 1;
    end
    coli_model{i_model}(k) = coli_data{i_model}(k);

    for m = k + 1:nData{i_model}
        % Determination of per mechanism decay rate
        k_dark_inst = k_dark*teta_dark^(temp_data{i_model}(m-1)-20);
        k_pH_inst = k_pH*teta_pH^(temp_data{i_model}(m-1)-20)*...
            10^(pH_data{i_model}(m-1)-14);
        k_sun_inst = alpha_sun*sun_data{i_model}(m-1)/...
            ((slope_sigma*TSSTable(i_model) + intercept_sigma)*d)*(1-exp(-(slope_sigma*TSSTable(i_model) + intercept_sigma)*d));
        % Determination of E. coli cell count
        coli_model{i_model}(m) = coli_model{i_model}(m-1)*...
            (1 - (k_dark_inst + k_pH_inst + k_sun_inst)*...
            (time_data{i_model}(m)-time_data{i_model}(m-1))) + ...
            (CIN(i_model)- coli_model{i_model}(m-1))*QIN(i_model)/(d*S)*...
            (time_data{i_model}(m)-time_data{i_model}(m-1));
        
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
                k_pH_inst = k_pH*teta_pH^(temp_local(i_subdiv)-20)*...
                    10^(pH_local(i_subdiv)-14);
                k_sun_inst = alpha_sun*sun_local(i_subdiv)/...
                    ((slope_sigma*TSSTable(i_model) + intercept_sigma)*d)*(1-exp(-(slope_sigma*TSSTable(i_model) + intercept_sigma)*d));
                % Determination of E. coli cell count
                coli_model_local(i_subdiv) = coli_model_local(i_subdiv-1)*...
                    (1 - (k_dark_inst + k_pH_inst + k_sun_inst)*...
                    (time_local(i_subdiv)-time_local(i_subdiv-1))) + ...
                    (CIN(i_model)- coli_model_local(i_subdiv-1))*...
                    QIN(i_model)/(d*S)*(time_local(i_subdiv)-time_local(i_subdiv-1));
            end
        coli_model{i_model}(m) = coli_model_local(n_subdiv+1);
        end
    end

    index_coli = find(~isnan(coli_data{i_model}));
    coli_fit = [coli_fit ; coli_model{i_model}(index_coli)];
    coli_fit_no1st = [coli_fit_no1st;[NaN ; coli_model{i_model}(index_coli(2:end))]]; 

end

nDataColi = length(coli_fit_no1st);
coli_model_bestFit = coli_model;
coli_predicted_logcounts = log10(coli_fit_no1st);

% Calculation of fitness parameters
MRAE = 100*mean(...
    abs(coli_logcounts(isnan(coli_predicted_logcounts) == 0) - coli_predicted_logcounts(isnan(coli_predicted_logcounts) == 0))...
    ./coli_logcounts(isnan(coli_predicted_logcounts) == 0));
R_squared = 1 - nansum((coli_predicted_logcounts-coli_logcounts).^2)/...
    sum((coli_logcounts(isnan(coli_predicted_logcounts) == 0) - ...
    mean(coli_logcounts(isnan(coli_predicted_logcounts) == 0))).^2);

% Drawing Figure S7.1
figure(11), clf, hold on
plot(coli_predicted_logcounts,coli_logcounts,'ob','MarkerSize',6,'LineWidth',2)
plot([2,3,4,5,6,7],[2,3,4,5,6,7],'-k','LineWidth',2)
errorbar(coli_predicted_logcounts,coli_logcounts,...
        coli_logcounts_IC_minus - coli_logcounts,...
        coli_logcounts_IC_plus - coli_logcounts,...
        '.k','MarkerSize',0.1,'LineWidth',1)
ylabel({'Measured {\itE. coli} cell counts' ; '(log_1_0 MPN{\cdot}100 mL)'},'FontSize',fs,...
    'FontWeight','bold','FontName','Arial')
xlabel({'Predicted {\itE. coli} cell counts' ; '(log_1_0 MPN{\cdot}100 mL)'},'FontSize',fs,...
    'FontWeight','bold','FontName','Arial')
ax2 = gca;
ax2.FontSize = fs;
ax2.XLim = [2 7];
ax2.YLim = [2 7];
title(strcat('R^2 = ',{' '},num2str(R_squared,formatSpec),{' '},';',...
    {' '},'MRAE =',{' '},num2str(MRAE,formatSpec),{' '},'%'),...
    'FontSize',fs,'FontWeight','bold','FontName','Arial')


fig = gcf;
fig.Position = [150 150 500 450];

%% Part IV Monte Carlo on fitness test


% This part aims at performing Monte Carlo analysis to draw the range of
% uncertainty on the fit of E. coli cell count due to input parameters
% uncertainty. Breifly, all input data are varied randomly within their
% uncertainty range (Table 2 in the manuscript) following normal laws and
% the fitting routine is carried out on the modified input dataset. This is
% repeated N_MC times, thus generating a large number of fitted parameters
% data sets. Theses sets of fitted parameters are then run to predict E.
% coli cell count. The results show how input uncertainty leads to
% uncertainty in the model fitness and enables to show the degree of
% uncertainty of the model fit.

% CALCULCATION OPTIONS

% option_MC: option allowing to choose to do a MC simulation (= 1), or not
    % (= 0). In the second case, it will set all errors to 0 and do a single
    % calculation. A last option is added to load the dataset from a
    % previous simulation (= 2).
% option_graphic: option to show the graphical output of Matlab
% optimization function.
    
option_MC = 2;
option_save_MC = 0;

% Monte Carlo initialisation

if option_MC == 0
    N_MC = 1;
else
    if option_MC == 1
        N_MC = 2000;
    else
        if option_MC == 2
            load('./Daily profiles/inputData_fit_MonteCarlo.mat')
        end
    end
end

% The simulation is ran on the validation data base
% NB; This simulation must be performed having already loaded results from
% the validation of the model (i.e. previous section)
dp_name = {'12_10' , '17_11' , '10_02' , '16_03' };
TSSTable = [220 , 310 , 210 , 740 ];
CIN = [3.47 , 6.15 , 4.96 , 3.45 ]*10^6;
QIN = [0.0984 , 0.1056 , 0.1488 , 0.12 ];
nExp = 4;

d = 0.25;
S = 3.42;
slope_sigma = 0.159;
intercept_sigma = 11.769;

% Fitted parameters values: the model fitting is initialized with the best
% fit values obtain for this study.

fitted_par = [k_pH_bestFit0,...
    teta_pH_bestFit0,...
    k_dark_bestFit0,...
    teta_dark_bestFit0,...
    alpha_sun_bestFit0];

% Uncertainty definition
if option_MC == 1
    err_pH = 0.2; % pH unit
    err_temp = 0.2; % °C
    err_sun = 0.1; % (relative, *100%)
    err_coli_IN = 0.15; % (relative, *100)
    err_coli = 1; % 1, -1, or value to compute it as basevalue+err_coli*(maxIC - minIC)/2
    err_TSS = 0.08; % (relative, *100%)
    err_slope_sigma = 0.0441; % (L/mg/m)
    err_intercept_sigma = 13.1; % m-1
    err_depth = 0.025; % m
    err_S = 0.171; % m2
    err_Q_IN = 0.1; % relative error
else
    if option_MC == 0
        err_pH = 0;
        err_temp = 0;
        err_sun = 0;
        err_coli_IN = 0;
        err_coli = 0;
        err_TSS = 0;
        err_slope_sigma = 0;
        err_intercept_sigma = 0;
        err_depth = 0;
        err_S = 0;
        err_Q_IN = 0;
    end
end

time_data = cell(nExp,1);
pH_data = cell(nExp,1);
temp_data = cell(nExp,1);
sun_data = cell(nExp,1);
coli_data = cell(nExp,1);
coli_data_IC_minus = cell(nExp,1);
coli_data_IC_plus = cell(nExp,1);
nData = cell(nExp,1);

for i = 1:nExp
    A = xlsread('./Daily profiles/Daily_profile_transformedData.xlsx',dp_name{i});
    
    time_data{i} = A(:,1);
    sun_data{i} = A(:,5);
    pH_data{i} = A(:,2);
    temp_data{i} = A(:,3);
    coli_data{i} = A(:,6);
    coli_data_IC_minus{i} = A(:,7);
    coli_data_IC_plus{i} = A(:,8);
    nData{i} = length(time_data{i});
end

% Monte Carlo loop

if option_MC ~= 2
    SSR_MC = NaN(1,N_MC);
    MRAE_MC = NaN(1,N_MC);
    coli_model_MC = cell(nExp,1); % Table of E. coli cell count predicted (MPN/100mL)
    for i = 1:nExp
        coli_model_MC{i} = NaN(nData{i},N_MC);
    end
    % Monte Carlo loop

    for i_MC = 1:N_MC
        % Definition of the fitting database
        [time_data_MC,sun_data_MC,pH_data_MC,temp_data_MC,coli_data_MC,TSS_MC,slope_sigma_MC,...
            intercept_sigma_MC,d_MC,S_MC,Q_IN_MC,C_IN_MC] =...
            create_database_MonteCarlo(time_data,sun_data,pH_data,temp_data,coli_data,...
                coli_data_IC_minus,coli_data_IC_plus,nData,...
                TSSTable,CIN,QIN,slope_sigma,intercept_sigma,d,S,...
                err_pH,err_temp,err_sun,err_coli_IN,err_coli,err_TSS,...
                err_slope_sigma,err_intercept_sigma,err_depth,err_S,err_Q_IN,...
                nExp);


        [SSR_MC(i_MC),MRAE_MC(i_MC),coli_model_loc,~] =...
            dailyProfile_fitting_function(fitted_par,...
            time_data_MC,sun_data_MC,pH_data_MC,temp_data_MC,coli_data_MC,nData,...
            TSS_MC,Q_IN_MC,C_IN_MC,slope_sigma_MC,intercept_sigma_MC,d_MC,S_MC,...
            nExp);
        for iExp = 1:nExp
            coli_model_MC{iExp,1}(:,i_MC) = coli_model_loc{iExp,1};
        end

    i_MC
    end
end

if option_save_MC == 1
    save('./Daily profiles/inputData_fit_MonteCarlo.mat',...
        'coli_data_MC','coli_model_MC','MRAE_MC','SSR_MC',...
        'pH_data_MC','sun_data_MC','temp_data_MC','time_data_MC',...
        'C_IN_MC','d_MC','intercept_sigma_MC','slope_sigma_MC','Q_IN_MC','S_MC','TSS_MC',...
        'N_MC')
end

% Uncertainty fit: creation of data tables for error plots

median_curve = cell(nExp,1);
area_table = cell(nExp,1);
prctile25_curve = cell(nExp,1);
prctile75_curve = cell(nExp,1);
prctile5_curve = cell(nExp,1);
prctile95_curve = cell(nExp,1);
max_curve = cell(nExp,1);
min_curve = cell(nExp,1);
for iExp = 1:nExp
    median_curve{iExp} = NaN(size(coli_model_MC{iExp},1),1);
    prctile25_curve{iExp} = NaN(size(coli_model_MC{iExp},1),1); 
    prctile75_curve{iExp} = NaN(size(coli_model_MC{iExp},1),1); 
    prctile5_curve{iExp} = NaN(size(coli_model_MC{iExp},1),1); 
    prctile95_curve{iExp} = NaN(size(coli_model_MC{iExp},1),1); 
    max_curve{iExp} = NaN(size(coli_model_MC{iExp},1),1);
    min_curve{iExp} = NaN(size(coli_model_MC{iExp},1),1);
    for iPoint = 1:size(coli_model_MC{iExp},1)
        median_curve{iExp}(iPoint) = prctile(log10(coli_model_MC{iExp,1}(iPoint,:)),50);
        prctile25_curve{iExp}(iPoint) = prctile(log10(coli_model_MC{iExp}(iPoint,:)),25);
        prctile75_curve{iExp}(iPoint) = prctile(log10(coli_model_MC{iExp}(iPoint,:)),75);
        prctile5_curve{iExp}(iPoint) = prctile(log10(coli_model_MC{iExp}(iPoint,:)),5);
        prctile95_curve{iExp}(iPoint) = prctile(log10(coli_model_MC{iExp}(iPoint,:)),95);
        max_curve{iExp}(iPoint) = max(log10(coli_model_MC{iExp}(iPoint,:)));
        min_curve{iExp}(iPoint) = min(log10(coli_model_MC{iExp}(iPoint,:)));
    end
    area_table{iExp} = [ min_curve{iExp} ,...
        prctile5_curve{iExp} - min_curve{iExp},...
        prctile25_curve{iExp} - prctile5_curve{iExp}, ...
        median_curve{iExp} - prctile25_curve{iExp},...
        prctile75_curve{iExp} - median_curve{iExp} ,...
        prctile95_curve{iExp} - prctile75_curve{iExp},...
        max_curve{iExp} - prctile95_curve{iExp}];
end

% Drawing of Figure 2

figure(30), clf, hold on
colours_MC = [1 1 1].*[1 ; 0.95 ; 0.75 ; 0.5 ; 0.5 ; 0.75; 0.95];
for i = 1:nExp
    subplot(2,2,i), hold on
    time_data_loc = datetime(0,1,floor(time_data{i}),0,0,3600*24*(time_data{i} -floor(time_data{i})));
    aplot = area(datenum(time_data_loc),area_table{i},1,'ShowBaseLine','off');
    for iColor = 1:7
        aplot(iColor).FaceColor = colours_MC(iColor,:);
        aplot(iColor).LineStyle = 'none';
    end
    bplot = plot(datenum(time_data_loc),log10(coli_model_bestFit{i}),'-k','LineWidth',2) ;
    cplot = plot(datenum(time_data_loc),log10(coli_data{i}),'ob','MarkerSize',6,'LineWidth',2) ;
    dplot = errorbar(datenum(time_data_loc),log10(coli_data{i}),...
        log10(coli_data{i}) - log10(coli_data_IC_minus{i}),...
        log10(coli_data_IC_plus{i}) - log10(coli_data{i}),...
        '.k','MarkerSize',0.1,'LineWidth',1);
%     eplot = plot(datenum(time_data_loc),log10(coli_model{i}),'--k','LineWidth',2);
    ylim([3 6])
    datetick('x')
    ax1 = gca;
    ax1.FontSize = fs-2;
    ax1.XLim = [datenum(dateshift(time_data_loc(1),'start','day') + 8/24),...
        datenum(dateshift(time_data_loc(1),'start','day') + 1 + 10/24)];
    xticks(ax1,[datenum(dateshift(time_data_loc(1),'start','day') + 9/24),...
        datenum(dateshift(time_data_loc(1),'start','day') + 21/24), ...
        datenum(dateshift(time_data_loc(1),'start','day') + 1 + 9/24)]);
    xticklabels(ax1,{datestr(dateshift(time_data_loc(1),'start','day') + 9/24,'dd/mm HH:MM'),...
        datestr(dateshift(time_data_loc(1),'start','day') + 21/24,'dd/mm HH:MM'), ...
        datestr(dateshift(time_data_loc(1),'start','day') + 1 + 9/24,'dd/mm HH:MM')});
    xtickangle(ax1,45)
    yticks(ax1,[3 4 5 6])
    set(ax1,'Box','off')
    if i <= 2
        title('Spring','FontSize',fs,'FontWeight','bold','FontName','Arial')
    else
        if i == 3
            title('Summer','FontSize',fs,'FontWeight','bold','FontName','Arial')
        else
            title('Fall','FontSize',fs,'FontWeight','bold','FontName','Arial')
        end
    end
    
    set(gca, 'Layer', 'top')
end
legend([aplot(2) aplot(3) aplot(4) bplot(1) cplot(1)],{'Outside 5/95 percentile','5-25 or 75-95 percentile',...
    '25 - 75 percentile','Best fit','Measured data'},'FontSize',fs-4,'FontWeight','bold',...
    'FontName','Arial','Location','SouthEast')
subplot(2,2,1)
ylabel({'{\itE. coli} cell counts' ; '(log_1_0 MPN{\cdot}100 mL^-^1)'},'FontSize',fs,...
    'FontWeight','bold','FontName','Arial')
subplot(2,2,3)
ylabel({'{\itE. coli} cell counts' ; '(log_1_0 MPN{\cdot}100 mL^-^1)'},'FontSize',fs,...
    'FontWeight','bold','FontName','Arial')

fig = gcf;
fig.Position = [50 50 750 700];
