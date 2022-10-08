%% Main_sensitivityAnalysis_fittedParameters

% This script aims at performing the fitting model sensitivity analysis and
% drawing the tornado diagram associated corresponding to Figure 1 of the
% manuscript "Modelling E. coli removal during real domestic wastewater
% treatment in outdoor pilot scale High Rate Algal Ponds".

% Briefly, each input parameters is alternatively varied to its max and min
% value of its uncertainty interval and the fitted parameters are computed
% as performed for the main calculation, except that the value of the best
% fit is used as initial value of the fitted parameters in the fitting
% routine.

%% Options
clear all
fs = 15;

% option_graphic: option to create a plot output assessing the quality of
% the fit (= 1)
option_graphic = 0;

%% Import of parameters

d = 0.25;
S = 3.42;
slope_sigma = 0.159;
intercept_sigma = 11.769;

n_param = 9;
N_SA = 2*n_param;


% The cal
nExp = 3;
dp_name = {'30_09' , '29_10' , '03_02'};
% light absorption by experiment (m-1)
TSSTable = [248 , 295 , 265 ];
% Inlet E. coli by experiment (MPN/100 mL)
CIN = [4.36 , 3.94 , 26 ]*10^6;
% Inlet flowrate by experiment (m3/d)
QIN = [0.1056 , 0.0984 , 0.144 ];

% Import of data

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

%% Best fit
% A first run to get the best fit values of the fitted parameters is
% performed.

% Fitting outputs definition

SSR_tot = NaN;
k_pH = NaN;
teta_pH = NaN;
k_dark = NaN;
teta_dark = NaN;
alpha_sun = NaN;
MRAE = NaN;

% Fit
    % The fitted parameters calculated in the best fit at bench scale are
    % used as the initial values of the fitting

k_pH_0 = 2860;
teta_pH_0 = 1.45;
k_dark_0 = 47.4;
teta_dark_0 = 1.00;
alpha_sun_0 = 0;

x0 = [k_pH_0,teta_pH_0,k_dark_0,teta_dark_0,alpha_sun_0];
lb = [0,1,0,1,0];
ub = [500000,10,200,10,10];
    
A = []; b = []; Aeq = []; beq = []; nonlcon = [];

options = optimoptions('fmincon','MaxIterations',200,'PlotFcn','optimplotfval');

[xfit,SSR,exitflag,output] = fmincon(@(x)dailyProfile_fitting_function_SSR(x,...
    time_data,sun_data,pH_data,temp_data,coli_data,nData,...
    TSSTable,QIN,CIN,slope_sigma,intercept_sigma,d,S,...
    nExp),...
    x0,...
    A,b,Aeq,beq,lb,ub,nonlcon,...
    options);

[SSR,MRAE_bestFit0] = dailyProfile_fitting_function(xfit,...
    time_data,sun_data,pH_data,temp_data,coli_data,nData,...
    TSSTable,QIN,CIN,slope_sigma,intercept_sigma,d,S,...
    nExp);

% New fitting parameters
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


%% Sensitivity analysis run
% SA Errors
% Errors below are given in 95% CI
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

SSR_tot_SA = NaN(1,N_SA);
MRAE_tot_SA = NaN(1,N_SA);
k_pH_SA = NaN(1,N_SA);
teta_pH_SA = NaN(1,N_SA);
k_dark_SA = NaN(1,N_SA);
teta_dark_SA = NaN(1,N_SA);
alpha_sun_SA = NaN(1,N_SA);

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
   
    k_pH_0 = k_pH_bestFit0;
    teta_pH_0 = teta_pH_bestFit0;
    k_dark_0 = k_dark_bestFit0;
    teta_dark_0 = teta_dark_bestFit0;
    alpha_sun_0 = alpha_sun_bestFit0;

% Calculation
    x0 = [k_pH_0,teta_pH_0,k_dark_0,teta_dark_0,alpha_sun_0];
    lb = [0,1,0,1,0];
    ub = [500000,10,200,10,10];
    
    A = []; b = []; Aeq = []; beq = []; nonlcon = [];
    if option_graphic == 1
        options = optimoptions('fmincon','MaxIterations',400,'PlotFcn','optimplotfval');
        options.MaxFunctionEvaluations = 5000;
    else
        options = optimoptions('fmincon','MaxIterations',400);
    end
    

    [xfit,SSR,exitflag,output] = fmincon(@(x)dailyProfile_fitting_function_SSR(x,...
        time_data_SA,sun_data_SA,pH_data_SA,temp_data_SA,coli_data_SA,nData,...
        TSS_SA,Q_IN_SA,C_IN_SA,slope_sigma_SA,intercept_sigma_SA,d_SA,S_SA,...
        nExp),...
        x0,...
        A,b,Aeq,beq,lb,ub,nonlcon,...
        options)

    [SSR,MRAE,coli_model_SA,log_coli_fit_SA] = dailyProfile_fitting_function(xfit,...
        time_data_SA,sun_data_SA,pH_data_SA,temp_data_SA,coli_data_SA,nData,...
        TSS_SA,Q_IN_SA,C_IN_SA,slope_sigma_SA,intercept_sigma_SA,d_SA,S_SA,...
        nExp);
   
    SSR_tot_SA(i_SA) = SSR;
    MRAE_tot_SA(i_SA) = MRAE;
    k_pH_SA(i_SA) = xfit(1);
    teta_pH_SA(i_SA) = xfit(2);
    k_dark_SA(i_SA) = xfit(3);
    teta_dark_SA(i_SA) = xfit(4);
    alpha_sun_SA(i_SA) = xfit(5);
  
end

%% Plot drawing

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

k_pH_high = k_pH_SA(1:2:end)' ;
k_pH_low = k_pH_SA(2:2:end)' ;

teta_pH_high =  teta_pH_SA(1:2:end)' ;
teta_pH_low = teta_pH_SA(2:2:end)' ;

k_dark_high = k_dark_SA(1:2:end)' ;
k_dark_low = k_dark_SA(2:2:end)' ;

teta_dark_high =  teta_dark_SA(1:2:end)' ;
teta_dark_low = teta_dark_SA(2:2:end)' ;

alpha_sun_high =  alpha_sun_SA(1:2:end)' ;
alpha_sun_low = alpha_sun_SA(2:2:end)' ;

MRAE_high = MRAE_tot_SA(1:2:end)';
MRAE_low = MRAE_tot_SA(2:2:end)';

% Base value

k_pH_BV = k_pH_bestFit0;
teta_pH_BV = teta_pH_bestFit0;
k_dark_BV = k_dark_bestFit0;
teta_dark_BV = teta_dark_bestFit0;
alpha_sun_BV = alpha_sun_bestFit0;
MRAE_BV = MRAE_bestFit0;


% Ordering of the tornado based on dark decay

B = [k_dark_high,k_dark_low];
C = abs(B - k_dark_BV);
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

k_pH_high = k_pH_high(index_swap,:);
k_pH_low = k_pH_low(index_swap,:);

teta_pH_high = teta_pH_high(index_swap,:);
teta_pH_low = teta_pH_low(index_swap,:);

k_dark_high = k_dark_high(index_swap,:);
k_dark_low = k_dark_low(index_swap,:);

teta_dark_high = teta_dark_high(index_swap,:);
teta_dark_low = teta_dark_low(index_swap,:);

alpha_sun_high =  alpha_sun_high(index_swap,:);
alpha_sun_low = alpha_sun_low(index_swap,:);

MRAE_high = MRAE_high(index_swap,:);
MRAE_low = MRAE_low(index_swap,:);

%% Plot of manuscript Figure 1

step_fig = 4;
X = 1:step_fig:step_fig*(n-1) + 2  ;


figure (2), clf, hold on

subplot(2,3,3), hold on
barh(X,(alpha_sun_high),'FaceColor',[0.75,0.75,0.75],'FaceAlpha',0.5,'BaseValue',alpha_sun_BV,'LineWidth',2)
barh(X,(alpha_sun_low),'FaceColor',[1,1,1],'FaceAlpha',0.5,'BaseValue',alpha_sun_BV,'LineWidth',2)
ax = gca; 
set(ax,'ytick',X,'yticklabel',[],'FontSize',fs,'FontWeight','bold','Position',[0.77 0.60 0.21 0.39])
ylim([0 step_fig*(n) + 1]);
xlabel('\alpha (m^2{\cdot}W^{-1}{\cdot}d^{-1})','FontSize',fs,'FontWeight','bold');


subplot(2,3,2), hold on
barh(X,(k_pH_high),'FaceColor',[0.75,0.75,0.75],'FaceAlpha',0.5,'BaseValue',k_pH_BV,'LineWidth',2)
barh(X,(k_pH_low),'FaceColor',[1,1,1],'FaceAlpha',0.5,'BaseValue',k_pH_BV,'LineWidth',2)
ax = gca; 
set(ax,'yticklabel',[],'FontSize',fs,'FontWeight','bold','Position',[0.54 0.60 0.21 0.39])
ylim([0 step_fig*(n) + 1]);
xlabel('k_{20}^{pH} ({\cdot}d^{-1})','FontSize',fs,'FontWeight','bold');
xtickformat('%.3f')


subplot(2,3,5), hold on
barh(X,(teta_pH_high),'FaceColor',[0.75,0.75,0.75],'FaceAlpha',0.5,'BaseValue',teta_pH_BV,'LineWidth',2)
barh(X,(teta_pH_low),'FaceColor',[1,1,1],'FaceAlpha',0.5,'BaseValue',teta_pH_BV,'LineWidth',2)
ax = gca; 
set(ax,'ytick',X,'yticklabel',[],'FontSize',fs,'FontWeight','bold','Position',[0.54, 0.09 0.21 0.39]);
ylim([0 step_fig*(n) + 1]);
xlabel('\theta^{pH}','FontSize',fs,'FontWeight','bold');



subplot(2,3,1), hold on
barh(X,(k_dark_high),'FaceColor',[0.75,0.75,0.75],'FaceAlpha',0.5,'BaseValue',k_dark_BV,'LineWidth',2)
barh(X,(k_dark_low),'FaceColor',[1,1,1],'FaceAlpha',0.5,'BaseValue',k_dark_BV,'LineWidth',2)
ax = gca; 
set(ax,'ytick',X,'yticklabel',Legend,'FontSize',fs,'FontWeight','bold','Position',[0.30 0.60 0.21 0.39])
ylim([0 step_fig*(n) + 1]);
xlabel('k_{20}^{dark} ({\cdot}d^{-1})','FontSize',fs+2,'FontWeight','bold');

subplot(2,3,4), hold on
barh(X,(teta_dark_high),'FaceColor',[0.75,0.75,0.75],'FaceAlpha',0.5,'BaseValue',teta_dark_BV,'LineWidth',2)
barh(X,(teta_dark_low),'FaceColor',[1,1,1],'FaceAlpha',0.5,'BaseValue',teta_dark_BV,'LineWidth',2)
ax = gca; 
set(ax,'ytick',X,'yticklabel',Legend,'FontSize',fs,'FontWeight','bold','Position',[0.30 0.09 0.21 0.39]); 
ylim([0 step_fig*(n) + 1]);
xlabel('\theta^{dark}','FontSize',fs,'FontWeight','bold');

fig = gcf;
fig.Position = [50 50 1110 700];

