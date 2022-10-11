%% Main_fittingParameters_MCanalysis

% This script aims at performing Monte Carlo analysis to determine a set of
% fitted parameters set in the publication 'Modelling E. coli removal during real
% domestic wastewater treatment in outdoor pilot scale High Rate Algal
% Ponds'. Each set of fitted paramaeters hence determined represent
% plausible values due to input uncertainty where the fitting was
% performed. 

% Briefly, values within uncertainty range of all parameters are generated
% and the fitting method as performed in Main_dailyProfile_analysis.m is
% repeated here each time. N_MC such fitting are performed here.


%% Options
clear all
fs = 14;
formatSpec = '%.2f'; % option to set only two digits in results used in figure titles


%% Part II: Fit of data based on only 3 daily profiles: 

dp_name = {'30_09' , '29_10' , '03_02'};
% light absorption by experiment (m-1)
TSSTable = [248 , 295 , 265 ];
% Inlet E. coli by experiment (MPN/100 mL)
CIN = [4.36 , 3.94 , 26 ]*10^6;
% Inlet flowrate by experiment (m3/d)
QIN = [0.1056 , 0.0984 , 0.144 ];
% depth (m)
d = 0.25;
% Volume (m3)
S = 3.42;
nExp = 3;

% Light absorption coefficient
slope_sigma = 0.159;
intercept_sigma = 11.769;


time_data = cell(nExp,1);
pH_data = cell(nExp,1);
temp_data = cell(nExp,1);
sun_data = cell(nExp,1);
coli_data = cell(nExp,1);
nData = cell(nExp,1);
coli_data_IC_minus = cell(nExp,1);
coli_data_IC_plus = cell(nExp,1);

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


% Fitting outputs definition

SSR_tot = NaN;
k_pH = NaN;
teta_pH = NaN;
k_dark = NaN;
teta_dark = NaN;
alpha_sun = NaN;


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

SSR_tot = SSR;

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


%% Uncertainty analysis - Monte Carlo

% This part aims at performing Monte Carlo analysis to determine the
% uncertainty associated to fitted parameters for the validated model as
% used in paper 3b based on the strategy 4 developed (i.e. using Matlab
% built-in function fminsearch).
 
% OPTIONS
 
option_save_MC = 0;
 
% Monte Carlo uncertainty initialisation
 
N_MC = 2000;

dp_name = {'30_09' , '29_10' , '03_02'};
% light absorption by experiment (m-1)
TSSTable = [248 , 295 , 265 ];
% Inlet E. coli by experiment (MPN/100 mL)
CIN = [4.36 , 3.94 , 26 ]*10^6;
% Inlet flowrate by experiment (m3/d)
QIN = [0.1056 , 0.0984 , 0.144 ];
% depth (m)
nExp = 3;
% depth (m)
d = 0.25;
% Surface (m2)
S = 3.42;

slope_sigma = 0.159;
intercept_sigma = 11.769;


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

% Monte Carlo outputs definition

SSR_tot_MC = NaN(1,N_MC);
MRAE_MC = NaN(1,N_MC);
k_pH_MC = NaN(1,N_MC);
teta_pH_MC = NaN(1,N_MC);
k_dark_MC = NaN(1,N_MC);
teta_dark_MC = NaN(1,N_MC);
alpha_sun_MC = NaN(1,N_MC);

% Monte Carlo loop

for i_MC = 1:N_MC
  %  Definition of the fitting database
    [time_data_MC,sun_data_MC,pH_data_MC,temp_data_MC,coli_data_MC,TSS_MC,slope_sigma_MC,...
        intercept_sigma_MC,d_MC,S_MC,Q_IN_MC,C_IN_MC] =...
        create_database_MonteCarlo(time_data,sun_data,pH_data,temp_data,coli_data,...
            coli_data_IC_minus,coli_data_IC_plus,nData,...
            TSSTable,CIN,QIN,slope_sigma,intercept_sigma,d,S,...
            err_pH,err_temp,err_sun,err_coli_IN,err_coli,err_TSS,...
            err_slope_sigma,err_intercept_sigma,err_depth,err_S,err_Q_IN,...
            nExp);

SSR_tot = SSR;
    % The fitted parameters calculated in the base simulation are used as the
    % initial values of the fitting
    
    
    
    
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
    options = optimoptions('fmincon','MaxIterations',100);
    
    [xfit,SSR,exitflag,output] = fmincon(@(x)dailyProfile_fitting_function_SSR(x,...
    time_data_MC,sun_data_MC,pH_data_MC,temp_data_MC,coli_data_MC,nData,...
    TSS_MC,Q_IN_MC,C_IN_MC,slope_sigma_MC,intercept_sigma_MC,d_MC,S_MC,...
    nExp),...
    x0,...
    A,b,Aeq,beq,lb,ub,nonlcon,...
    options);


    SSR_tot_MC(i_MC) = SSR;
    k_pH_MC(i_MC) = xfit(1);
    teta_pH_MC(i_MC) = xfit(2);
    k_dark_MC(i_MC) = xfit(3);
    teta_dark_MC(i_MC) = xfit(4);
    alpha_sun_MC(i_MC) = xfit(5);
    
    i_MC/N_MC*100

end

if option_save_MC == 1
    save('./Daily profiles/fittingParametersMonteCarlo.mat','SSR_tot_MC','k_pH_MC','teta_pH_MC','k_dark_MC','teta_dark_MC','alpha_sun_MC')
end
