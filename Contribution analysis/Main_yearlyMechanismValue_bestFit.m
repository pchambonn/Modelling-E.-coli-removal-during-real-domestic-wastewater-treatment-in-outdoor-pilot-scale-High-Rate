%% Main_yearlyMechanismValue_bestFit

% THis script aims at computing the decay VALUE of each mechanism to
% overall E. coli decay in a theoretical HRAP which environmental
% conditions were as measured i HRAP A and B. These vectors serve as the
% basis for plotting Figure 4 and Figure S10.1 of the publication
% 'Modelling E. coli removal during real domestic wastewater treatment in
% outdoor pilot scale High Rate Algal Ponds' (as performed in the script
% Main_plot_seasonalZoomedDecay). 

% The contribution array store the ratio between specific mechnaism decay
% rate (d-1) over total sum of the decay rates. It is a 5x1 cell, each entry containing a the data set of values of a season:
    % Entry 1 is for winter,
    % Entry 2 is for spring,
    % Entry 3 is for summer,
    % Entry 4 is for fall.
    % Entry 5 is for the total. NB: to be thouhgt about due to imbalance
        % between seasonal data set size, I may need to correct the data?
    
clear all

%%
load('./Daily profiles/fittingParameters_bestFit0.mat')
load('./Yearly environmental data/yearly_environmental_data_formatted.mat')
k_pH = k_pH_bestFit0;
teta_pH = teta_pH_bestFit0;
k_dark = k_dark_bestFit0;
teta_dark = teta_dark_bestFit0;
alpha_sun = alpha_sun_bestFit0;

%%
n_A_winter = size(data_A_winter,1);
n_A_spring = size(data_A_spring,1);
n_A_summer = size(data_A_summer,1);
n_A_fall = size(data_A_fall,1);
n_A_total = n_A_winter + n_A_spring + n_A_summer + n_A_fall;

n_B_winter = size(data_B_winter,1);
n_B_spring = size(data_B_spring,1);
n_B_summer = size(data_B_summer,1);
n_B_fall = size(data_B_fall,1);
n_B_total = n_B_winter + n_B_spring + n_B_summer + n_B_fall;

d = 0.25;

% THe contribution array store the ratio between specific mechnaism decay
% rate (d-1) over total sum of the decay rates. 'Column' 1 is for winter, 2
% for spring, 3 for summer, and 4 for fall.

decay_pH_A = cell(5,1);
    decay_pH_A{1} = NaN(n_A_winter,1);
    decay_pH_A{2} = NaN(n_A_spring,1);
    decay_pH_A{3} = NaN(n_A_summer,1);
    decay_pH_A{4} = NaN(n_A_fall,1);
    decay_pH_A{5} = NaN(n_A_total,1);

decay_dark_A = cell(5,1);
    decay_dark_A{1} = NaN(n_A_winter,1);
    decay_dark_A{2} = NaN(n_A_spring,1);
    decay_dark_A{3} = NaN(n_A_summer,1);
    decay_dark_A{4} = NaN(n_A_fall,1);
    decay_dark_A{5} = NaN(n_A_total,1);

decay_sun_A = cell(5,1);
    decay_sun_A{1} = NaN(n_A_winter,1);
    decay_sun_A{2} = NaN(n_A_spring,1);
    decay_sun_A{3} = NaN(n_A_summer,1);
    decay_sun_A{4} = NaN(n_A_fall,1);
    decay_sun_A{5} = NaN(n_A_total,1);

decay_pH_B = cell(5,1);
    decay_pH_B{1} = NaN(n_B_winter,1);
    decay_pH_B{2} = NaN(n_B_spring,1);
    decay_pH_B{3} = NaN(n_B_summer,1);
    decay_pH_B{4} = NaN(n_B_fall,1);
    decay_pH_B{5} = NaN(n_B_total,1);

decay_dark_B = cell(5,1);
    decay_dark_B{1} = NaN(n_B_winter,1);
    decay_dark_B{2} = NaN(n_B_spring,1);
    decay_dark_B{3} = NaN(n_B_summer,1);
    decay_dark_B{4} = NaN(n_B_fall,1);
    decay_dark_B{5} = NaN(n_B_total,1);
        
decay_sun_B = cell(5,1);
    decay_sun_B{1} = NaN(n_B_winter,1);
    decay_sun_B{2} = NaN(n_B_spring,1);
    decay_sun_B{3} = NaN(n_B_summer,1);
    decay_sun_B{4} = NaN(n_B_fall,1);
    decay_sun_B{5} = NaN(n_B_total,1);

    
%% 
    % Simulation HRAP A

for i = 1:n_A_winter
    decay_pH_A{1}(i) = k_pH*10^(data_A_winter(i,2)-14)*teta_pH^(data_A_winter(i,3)-20);
    decay_dark_A{1}(i) = k_dark*teta_dark^(data_A_winter(i,3)-20);
    decay_sun_A{1}(i) = alpha_sun*data_A_winter(i,5)/(data_A_winter(i,4)*d)*(1 - exp(-data_A_winter(i,4)*d));
 end

for i = 1:n_A_spring
    decay_pH_A{2}(i) = k_pH*10^(data_A_spring(i,2)-14)*teta_pH^(data_A_spring(i,3)-20);
    decay_dark_A{2}(i) = k_dark*teta_dark^(data_A_spring(i,3)-20);
    decay_sun_A{2}(i) = alpha_sun*data_A_spring(i,5)/(data_A_spring(i,4)*d)*(1 - exp(-data_A_spring(i,4)*d));
end

for i = 1:n_A_summer
    decay_pH_A{3}(i) = k_pH*10^(data_A_summer(i,2)-14)*teta_pH^(data_A_summer(i,3)-20);
    decay_dark_A{3}(i) = k_dark*teta_dark^(data_A_summer(i,3)-20);
    decay_sun_A{3}(i) = alpha_sun*data_A_summer(i,5)/(data_A_summer(i,4)*d)*(1 - exp(-data_A_summer(i,4)*d));
end

for i = 1:n_A_fall       
    decay_pH_A{4}(i) = k_pH*10^(data_A_fall(i,2)-14)*teta_pH^(data_A_fall(i,3)-20);
    decay_dark_A{4}(i) =  k_dark*teta_dark^(data_A_fall(i,3)-20);
    decay_sun_A{4}(i) = alpha_sun*data_A_fall(i,5)/(data_A_fall(i,4)*d)*(1 - exp(-data_A_fall(i,4)*d));
end
    
decay_pH_A{5} = [decay_pH_A{1} ; decay_pH_A{2} ; decay_pH_A{3}; decay_pH_A{4}];
decay_dark_A{5} = [decay_dark_A{1} ; decay_dark_A{2} ; decay_dark_A{3}; decay_dark_A{4}];
decay_sun_A{5} = [decay_sun_A{1} ; decay_sun_A{2} ; decay_sun_A{3}; decay_sun_A{4}];


% Calculation HRAP B

 for i = 1:n_B_winter
    decay_pH_B{1}(i) = k_pH*10^(data_B_winter(i,2)-14)*teta_pH^(data_B_winter(i,3)-20);
    decay_dark_B{1}(i) = k_dark*teta_dark^(data_B_winter(i,3)-20);
    decay_sun_B{1}(i) = alpha_sun*data_B_winter(i,5)/(data_B_winter(i,4)*d)*(1 - exp(-data_B_winter(i,4)*d));
 end

for i = 1:n_B_spring
    decay_pH_B{2}(i) = k_pH*10^(data_B_spring(i,2)-14)*teta_pH^(data_B_spring(i,3)-20);
    decay_dark_B{2}(i) = k_dark*teta_dark^(data_B_spring(i,3)-20);
    decay_sun_B{2}(i) = alpha_sun*data_B_spring(i,5)/(data_B_spring(i,4)*d)*(1 - exp(-data_B_spring(i,4)*d));
end

for i = 1:n_B_summer
    decay_pH_B{3}(i) = k_pH*10^(data_B_summer(i,2)-14)*teta_pH^(data_B_summer(i,3)-20);
    decay_dark_B{3}(i) = k_dark*teta_dark^(data_B_summer(i,3)-20);
    decay_sun_B{3}(i) = alpha_sun*data_B_summer(i,5)/(data_B_summer(i,4)*d)*(1 - exp(-data_B_summer(i,4)*d));
end

for i = 1:n_B_fall       
    decay_pH_B{4}(i) = k_pH*10^(data_B_fall(i,2)-14)*teta_pH^(data_B_fall(i,3)-20);
    decay_dark_B{4}(i) =  k_dark*teta_dark^(data_B_fall(i,3)-20);
    decay_sun_B{4}(i) = alpha_sun*data_B_fall(i,5)/(data_B_fall(i,4)*d)*(1 - exp(-data_B_fall(i,4)*d));
end

decay_pH_B{5} = [decay_pH_B{1} ; decay_pH_B{2} ; decay_pH_B{3}; decay_pH_B{4}];
decay_dark_B{5} = [decay_dark_B{1} ; decay_dark_B{2} ; decay_dark_B{3}; decay_dark_B{4}];
decay_sun_B{5} = [decay_sun_B{1} ; decay_sun_B{2} ; decay_sun_B{3}; decay_sun_B{4}];


%%

save(strcat('./Contribution analysis/yearlyDecayValue.mat'),...
    'decay_pH_A','decay_dark_A','decay_sun_A',...
    'decay_pH_B','decay_dark_B','decay_sun_B');
