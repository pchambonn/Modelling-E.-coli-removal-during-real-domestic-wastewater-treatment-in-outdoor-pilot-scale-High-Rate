%% Main_sensitivityAnalysis_model_fitness

% This script aims at performing the sensitivity analysis of model fitness
% and drawing the tornado diagram associated (Figure S8-1) for  the
% manuscript "Modelling E. coli removal during real domestic wastewater
% treatment in outdoor pilot scale High Rate Algal Ponds".

% Briefly, each input parameters is alternatively varied to its max and min
% value of its uncertainty interval and the fitness outputs (MRAE and SSR)
% are computed for each varied input dataset.

%% Options
clear all
fs = 15;

% option_graphic: option to create outputs assessing the quality of the fit
    % (= 1)
option_graphic = 0;


%% Import parameters
n_param = 9;
N_SA = 2*n_param;


d = 0.25;
S = 3.42;
slope_sigma = 0.159;
intercept_sigma = 11.769;
dp_name = {'12_10' , '17_11' , '10_02' , '16_03' };
% TSS by experiment (mg/L)
TSSTable = [220 , 310 , 210 , 740 ];
% Inlet E. coli by experiment (MPN/100 mL)
CIN = [3.47 , 6.15 , 4.96 , 3.45 ]*10^6;
% Inlet flowrate by experiment (m3/d)
QIN = [0.0984 , 0.1056 , 0.1488 , 0.12 ];
nExp = 4;


% Import data

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

% Fitted parameters values
fitted_par = [3.5317e+04,...
    1.5531,...
    13.5605,...
    1.1139,...
    0.1314];
%% Base value from best fit

% Fitting outputs definition

[SSR_bestFit0,MRAE_bestFit0] = dailyProfile_fitting_function(fitted_par,...
    time_data,sun_data,pH_data,temp_data,coli_data,nData,...
    TSSTable,QIN,CIN,slope_sigma,intercept_sigma,d,S,...
    nExp);


%% Sensitivity analysis run
% SA Errors
% Errors below are given in 95% CI (relative or absolute value)

err_pH = 0.2; % pH unit
err_temp = 0.2; % °C
err_sun = 0.1; % (relative, *100%)
err_coli_IN = 0.15; % (relative, *100)
err_coli = 0; % 1, -1, or value to compute it as basevalue+err_coli*(maxIC - minIC)/2
err_TSS = 0.08; % (relative, *100%)
err_slope_sigma = 0.0441; % (L/mg/m)
err_intercept_sigma = 13.1; % m-1
err_depth = 0.025; % m
err_S = 0.171; % m2
err_Q_IN = 0.1; % relative error

error_matrix = zeros(n_param,N_SA);
for i = 1:n_param
    error_matrix(i,2*i-1) = 1;
    error_matrix(i,2*i) = -1;
end

% SA outputs definition

SSR_SA = NaN(1,N_SA);
MRAE_SA = NaN(1,N_SA);

% SA calculation loop

for i_SA = 1:N_SA
    % Definition of the fitting database
    [time_data_SA,sun_data_SA,pH_data_SA,temp_data_SA,coli_data_SA,TSS_SA,slope_sigma_SA,intercept_sigma_SA,d_SA,S_SA,Q_IN_SA,C_IN_SA] =...
        create_database_SensitivityAnalysis(time_data,sun_data,pH_data,temp_data,coli_data,coli_data_IC_minus,coli_data_IC_plus,nData,...
        TSSTable,CIN,QIN,slope_sigma,intercept_sigma,d,S,...
        err_pH*error_matrix(1,i_SA),err_temp*error_matrix(2,i_SA),err_sun*error_matrix(3,i_SA),...
        err_coli_IN*error_matrix(4,i_SA),err_coli,err_TSS*error_matrix(5,i_SA),...
        err_slope_sigma*error_matrix(6,i_SA),err_intercept_sigma*error_matrix(6,i_SA),... % sensitivity to slope and intercept are computed together because only trying over or underestimation of the slope.
        err_depth*error_matrix(7,i_SA),err_S*error_matrix(8,i_SA),err_Q_IN*error_matrix(9,i_SA),...
        nExp);
   
    
    [SSR,MRAE,coli_model_SA,log_coli_fit_SA] = dailyProfile_fitting_function(fitted_par,...
        time_data_SA,sun_data_SA,pH_data_SA,temp_data_SA,coli_data_SA,nData,...
        TSS_SA,Q_IN_SA,C_IN_SA,slope_sigma_SA,intercept_sigma_SA,d_SA,S_SA,...
        nExp);
   
    SSR_SA(i_SA) = SSR;
    MRAE_SA(i_SA) = MRAE;
    
end

%% Drawing of Figure S8.1

% This plot aims at drawing the tornado diagram for paper 4 model
% sensitivity analysis on the determination of fitted parameters

Y1 = 'Broth pH' ;
Y2 = 'Broth temperature';
Y3 = 'Sunlight intensity';
Y4 = '{\itE. coli} cell count in the influent';
Y5 = 'TSS';
Y6 = 'Light absorption coefficient';
Y7 = 'Reactor depth';
Y8 = 'Reactor surface';
Y9 = 'Influent flowrate';

Legend_xls = {Y1;Y2;Y3;Y4;Y5;Y6;Y7;Y8;Y9};

% Vectors of data

MRAE_high = MRAE_SA(1:2:end)' ;
MRAE_low = MRAE_SA(2:2:end)' ;

SSR_high = SSR_SA(1:2:end)';
SSR_low = SSR_SA(2:2:end)';

% Base value

SSR_BV = SSR_bestFit0;
MRAE_BV = MRAE_bestFit0;

% Ordering of the tornado

% Dark decay

B = [MRAE_high,MRAE_low];
C = abs(B - MRAE_BV);
n = size(C,1);
C_bis = [zeros(n,1),C];
for i = 1:n
    C_bis(i,1) = max(C_bis(i,2),C_bis(i,3));
end
[C_bis,index_swap] = sortrows(C_bis,1);

n = length(index_swap);
Legend = cell(n,1);
for i = 1:n
    Legend{i} = Legend_xls{index_swap(i)};
end


MRAE_high = MRAE_high(index_swap,:);
MRAE_low = MRAE_low(index_swap,:);

SSR_high = SSR_high(index_swap,:);
SSR_low = SSR_low(index_swap,:);


step_fig = 4;
X = 1:step_fig:step_fig*(n-1) + 2  ;

figure (2), clf, hold on

subplot(1,2,2), hold on
barh(X,(MRAE_high)*100,'FaceColor',[0.75,0.75,0.75],'FaceAlpha',0.5,'BaseValue',MRAE_BV*100,'LineWidth',2)
barh(X,(MRAE_low)*100,'FaceColor',[1,1,1],'FaceAlpha',0.5,'BaseValue',MRAE_BV*100,'LineWidth',2)
ax = gca; 
set(ax,'ytick',X,'yticklabel',[],'FontSize',fs,'FontWeight','bold','Position',[0.65 0.1 0.30 0.88])
ylim([0 step_fig*(n) + 1]);
xlabel('MRAE (%)','FontSize',fs,'FontWeight','bold');

subplot(1,2,1), hold on
barh(X,(SSR_high),'FaceColor',[0.75,0.75,0.75],'FaceAlpha',0.5,'BaseValue',SSR_BV,'LineWidth',2)
barh(X,(SSR_low),'FaceColor',[1,1,1],'FaceAlpha',0.5,'BaseValue',SSR_BV,'LineWidth',2)
ax = gca; 
set(ax,'ytick',X,'yticklabel',Legend,'FontSize',fs,'FontWeight','bold','Position',[0.30 0.1 0.30 0.88])
ylim([0 step_fig*(n) + 1]);
xlabel('SSR (log_1_0 MPN^2{\cdot}100 mL^-^2)','FontSize',fs,'FontWeight','bold');


fig = gcf;
fig.Position = [50 50 1020 720];

