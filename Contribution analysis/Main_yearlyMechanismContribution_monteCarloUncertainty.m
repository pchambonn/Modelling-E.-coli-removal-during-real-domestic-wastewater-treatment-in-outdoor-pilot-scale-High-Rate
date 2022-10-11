%% Main_yearlyMechanismContribution_monteCarloUncertainty

% THis script aims at computing the contribution of each mechanism to
% overall E. coli decay in a theoretical HRAP which environmental
% conditions were as measured i HRAP A and B including variations for every
% Monte Carlo set of value for fitting parameters. These vectors serve as
% the basis for plotting Figure S9.2 of the publication 'Modelling E. coli
% removal during real domestic wastewater treatment in outdoor pilot scale
% High Rate Algal Ponds' (as performed in the script
% Main_plot_yearlySeasonalContribution). 

% The contribution array store the ratio between specific mechnaism decay
% rate (d-1) over total sum of the decay rates. It is a 4x1 cell, each
% entry containing a the data set of values of a season:
    % Entry 1 is for winter,
    % Entry 2 is for spring,
    % Entry 3 is for summer,
    % Entry 4 is for fall.
    % Entry 5 is for the total.
% Each of these entries contains a N_data_season x N_MC double
% corresponding to the value variations for each Monte Carlo set of fitting
% parameters (per column) for each time stamp (per line).

clear all

option_save = 0;
%%
load('./Daily profiles/fittingParametersMonteCarlo.mat')
load('./Yearly environmental data/yearly_environmental_data_formatted.mat')

%% Variables needed to interpret from data imported

N_MC = length(k_pH_MC);

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


contribution_pH_A_MC = cell(5,1);
    contribution_pH_A_MC{1} = NaN(n_A_winter,N_MC);
    contribution_pH_A_MC{2} = NaN(n_A_spring,N_MC);
    contribution_pH_A_MC{3} = NaN(n_A_summer,N_MC);
    contribution_pH_A_MC{4} = NaN(n_A_fall,N_MC);
    contribution_pH_A_MC{5} = NaN(n_A_total,N_MC);


contribution_dark_A_MC = cell(5,1);
    contribution_dark_A_MC{1} = NaN(n_A_winter,N_MC);
    contribution_dark_A_MC{2} = NaN(n_A_spring,N_MC);
    contribution_dark_A_MC{3} = NaN(n_A_summer,N_MC);
    contribution_dark_A_MC{4} = NaN(n_A_fall,N_MC);
    contribution_dark_A_MC{5} = NaN(n_A_total,N_MC);

contribution_sun_A_MC = cell(5,1);
    contribution_sun_A_MC{1} = NaN(n_A_winter,N_MC);
    contribution_sun_A_MC{2} = NaN(n_A_spring,N_MC);
    contribution_sun_A_MC{3} = NaN(n_A_summer,N_MC);
    contribution_sun_A_MC{4} = NaN(n_A_fall,N_MC);
    contribution_sun_A_MC{5} = NaN(n_A_total,N_MC);

contribution_pH_B_MC = cell(5,1);
    contribution_pH_B_MC{1} = NaN(n_B_winter,N_MC);
    contribution_pH_B_MC{2} = NaN(n_B_spring,N_MC);
    contribution_pH_B_MC{3} = NaN(n_B_summer,N_MC);
    contribution_pH_B_MC{4} = NaN(n_B_fall,N_MC);
    contribution_pH_B_MC{5} = NaN(n_B_total,N_MC);

contribution_dark_B_MC = cell(5,1);
    contribution_dark_B_MC{1} = NaN(n_B_winter,N_MC);
    contribution_dark_B_MC{2} = NaN(n_B_spring,N_MC);
    contribution_dark_B_MC{3} = NaN(n_B_summer,N_MC);
    contribution_dark_B_MC{4} = NaN(n_B_fall,N_MC);
    contribution_dark_B_MC{5} = NaN(n_B_total,N_MC);

contribution_sun_B_MC = cell(5,1);
    contribution_sun_B_MC{1} = NaN(n_B_winter,N_MC);
    contribution_sun_B_MC{2} = NaN(n_B_spring,N_MC);
    contribution_sun_B_MC{3} = NaN(n_B_summer,N_MC);
    contribution_sun_B_MC{4} = NaN(n_B_fall,N_MC);
    contribution_sun_B_MC{5} = NaN(n_B_total,N_MC);


%%      
for i_MC = 1:N_MC
    k_pH = k_pH_MC(i_MC);
    teta_pH = teta_pH_MC(i_MC);
    k_dark = k_dark_MC(i_MC);
    teta_dark = teta_dark_MC(i_MC);
    alpha_sun = alpha_sun_MC(i_MC);
    
    % Simulation HRAP A
    
    decay_pH_A_winter = NaN(n_A_winter,1);
    decay_dark_A_winter = NaN(n_A_winter,1);
    decay_sun_A_winter = NaN(n_A_winter,1);
    
    decay_pH_A_spring = NaN(n_A_spring,1);
    decay_dark_A_spring = NaN(n_A_spring,1);
    decay_sun_A_spring = NaN(n_A_spring,1);
    
    decay_pH_A_summer = NaN(n_A_summer,1);
    decay_dark_A_summer = NaN(n_A_summer,1);
    decay_sun_A_summer = NaN(n_A_summer,1);
    
    decay_pH_A_fall = NaN(n_A_fall,1);
    decay_dark_A_fall = NaN(n_A_fall,1);
    decay_sun_A_fall = NaN(n_A_fall,1);
    
    for i = 1:n_A_winter
        decay_pH_A_winter(i) = k_pH*10^(data_A_winter(i,2)-14)*teta_pH^(data_A_winter(i,3)-20);
        decay_dark_A_winter(i) = k_dark*teta_dark^(data_A_winter(i,3)-20);
        decay_sun_A_winter(i) = alpha_sun*data_A_winter(i,5)/(data_A_winter(i,4)*d)*(1 - exp(-data_A_winter(i,4)*d));
        
        contribution_pH_A_MC{1}(i,i_MC) = decay_pH_A_winter(i)./(decay_pH_A_winter(i) + decay_dark_A_winter(i) + decay_sun_A_winter(i));
        contribution_dark_A_MC{1}(i,i_MC) = decay_dark_A_winter(i)./(decay_pH_A_winter(i) + decay_dark_A_winter(i) + decay_sun_A_winter(i));
        contribution_sun_A_MC{1}(i,i_MC) = decay_sun_A_winter(i)./(decay_pH_A_winter(i) + decay_dark_A_winter(i) + decay_sun_A_winter(i));
    end
    
    for i = 1:n_A_spring
        decay_pH_A_spring(i) = k_pH*10^(data_A_spring(i,2)-14)*teta_pH^(data_A_spring(i,3)-20);
        decay_dark_A_spring(i) = k_dark*teta_dark^(data_A_spring(i,3)-20);
        decay_sun_A_spring(i) = alpha_sun*data_A_spring(i,5)/(data_A_spring(i,4)*d)*(1 - exp(-data_A_spring(i,4)*d));
        
        contribution_pH_A_MC{2}(i,i_MC) = decay_pH_A_spring(i)./(decay_pH_A_spring(i) + decay_dark_A_spring(i) + decay_sun_A_spring(i));
        contribution_dark_A_MC{2}(i,i_MC) = decay_dark_A_spring(i)./(decay_pH_A_spring(i) + decay_dark_A_spring(i) + decay_sun_A_spring(i));
        contribution_sun_A_MC{2}(i,i_MC) = decay_sun_A_spring(i)./(decay_pH_A_spring(i) + decay_dark_A_spring(i) + decay_sun_A_spring(i));
    end
    
    for i = 1:n_A_summer
        decay_pH_A_summer(i) = k_pH*10^(data_A_summer(i,2)-14)*teta_pH^(data_A_summer(i,3)-20);
        decay_dark_A_summer(i) = k_dark*teta_dark^(data_A_summer(i,3)-20);
        decay_sun_A_summer(i) = alpha_sun*data_A_summer(i,5)/(data_A_summer(i,4)*d)*(1 - exp(-data_A_summer(i,4)*d));
        
        contribution_pH_A_MC{3}(i,i_MC) = decay_pH_A_summer(i)./(decay_pH_A_summer(i) + decay_dark_A_summer(i) + decay_sun_A_summer(i));
        contribution_dark_A_MC{3}(i,i_MC) = decay_dark_A_summer(i)./(decay_pH_A_summer(i) + decay_dark_A_summer(i) + decay_sun_A_summer(i));
        contribution_sun_A_MC{3}(i,i_MC) = decay_sun_A_summer(i)./(decay_pH_A_summer(i) + decay_dark_A_summer(i) + decay_sun_A_summer(i));
    end
    
    for i = 1:n_A_fall
        decay_pH_A_fall(i) = k_pH*10^(data_A_fall(i,2)-14)*teta_pH^(data_A_fall(i,3)-20);
        decay_dark_A_fall(i) = k_dark*teta_dark^(data_A_fall(i,3)-20);
        decay_sun_A_fall(i) = alpha_sun*data_A_fall(i,5)/(data_A_fall(i,4)*d)*(1 - exp(-data_A_fall(i,4)*d));
        
        contribution_pH_A_MC{4}(i,i_MC) = decay_pH_A_fall(i)./(decay_pH_A_fall(i) + decay_dark_A_fall(i) + decay_sun_A_fall(i));
        contribution_dark_A_MC{4}(i,i_MC) = decay_dark_A_fall(i)./(decay_pH_A_fall(i) + decay_dark_A_fall(i) + decay_sun_A_fall(i));
        contribution_sun_A_MC{4}(i,i_MC) = decay_sun_A_fall(i)./(decay_pH_A_fall(i) + decay_dark_A_fall(i) + decay_sun_A_fall(i));        
    end
    
    contribution_pH_A_MC{5}(:,i_MC) = ...
        [contribution_pH_A_MC{1}(:,i_MC);...
        contribution_pH_A_MC{2}(:,i_MC);...
        contribution_pH_A_MC{3}(:,i_MC);...
        contribution_pH_A_MC{4}(:,i_MC)];
    contribution_dark_A_MC{5}(:,i_MC) = ...
        [contribution_dark_A_MC{1}(:,i_MC);...
        contribution_dark_A_MC{2}(:,i_MC);...
        contribution_dark_A_MC{3}(:,i_MC);...
        contribution_dark_A_MC{4}(:,i_MC)];
    contribution_sun_A_MC{5}(:,i_MC) =...
        [contribution_sun_A_MC{1}(:,i_MC);...
        contribution_sun_A_MC{2}(:,i_MC);...
        contribution_sun_A_MC{3}(:,i_MC);...
        contribution_sun_A_MC{4}(:,i_MC)];
    
    % Calculation HRAP B
    
    decay_pH_B_winter = NaN(n_B_winter,1);
    decay_dark_B_winter = NaN(n_B_winter,1);
    decay_sun_B_winter = NaN(n_B_winter,1);
    
    decay_pH_B_spring = NaN(n_B_spring,1);
    decay_dark_B_spring = NaN(n_B_spring,1);
    decay_sun_B_spring = NaN(n_B_spring,1);
    
    decay_pH_B_summer = NaN(n_B_summer,1);
    decay_dark_B_summer = NaN(n_B_summer,1);
    decay_sun_B_summer = NaN(n_B_summer,1);
    
    decay_pH_B_fall = NaN(n_B_fall,1);
    decay_dark_B_fall = NaN(n_B_fall,1);
    decay_sun_B_fall = NaN(n_B_fall,1);
    
    for i = 1:n_B_winter
        decay_pH_B_winter(i) = k_pH*10^(data_B_winter(i,2)-14)*teta_pH^(data_B_winter(i,3)-20);
        decay_dark_B_winter(i) = k_dark*teta_dark^(data_B_winter(i,3)-20);
        decay_sun_B_winter(i) = alpha_sun*data_B_winter(i,5)/(data_B_winter(i,4)*d)*(1 - exp(-data_B_winter(i,4)*d));
        
        contribution_pH_B_MC{1}(i,i_MC) = decay_pH_B_winter(i)./(decay_pH_B_winter(i) + decay_dark_B_winter(i) + decay_sun_B_winter(i));
        contribution_dark_B_MC{1}(i,i_MC) = decay_dark_B_winter(i)./(decay_pH_B_winter(i) + decay_dark_B_winter(i) + decay_sun_B_winter(i));
        contribution_sun_B_MC{1}(i,i_MC) = decay_sun_B_winter(i)./(decay_pH_B_winter(i) + decay_dark_B_winter(i) + decay_sun_B_winter(i));
    end
    
    for i = 1:n_B_spring
        decay_pH_B_spring(i) = k_pH*10^(data_B_spring(i,2)-14)*teta_pH^(data_B_spring(i,3)-20);
        decay_dark_B_spring(i) = k_dark*teta_dark^(data_B_spring(i,3)-20);
        decay_sun_B_spring(i) = alpha_sun*data_B_spring(i,5)/(data_B_spring(i,4)*d)*(1 - exp(-data_B_spring(i,4)*d));
        
        contribution_pH_B_MC{2}(i,i_MC) = decay_pH_B_spring(i)./(decay_pH_B_spring(i) + decay_dark_B_spring(i) + decay_sun_B_spring(i));
        contribution_dark_B_MC{2}(i,i_MC) = decay_dark_B_spring(i)./(decay_pH_B_spring(i) + decay_dark_B_spring(i) + decay_sun_B_spring(i));
        contribution_sun_B_MC{2}(i,i_MC) = decay_sun_B_spring(i)./(decay_pH_B_spring(i) + decay_dark_B_spring(i) + decay_sun_B_spring(i));
    end
    
    for i = 1:n_B_summer
        decay_pH_B_summer(i) = k_pH*10^(data_B_summer(i,2)-14)*teta_pH^(data_B_summer(i,3)-20);
        decay_dark_B_summer(i) = k_dark*teta_dark^(data_B_summer(i,3)-20);
        decay_sun_B_summer(i) = alpha_sun*data_B_summer(i,5)/(data_B_summer(i,4)*d)*(1 - exp(-data_B_summer(i,4)*d));
        
        contribution_pH_B_MC{3}(i,i_MC) = decay_pH_B_summer(i)./(decay_pH_B_summer(i) + decay_dark_B_summer(i) + decay_sun_B_summer(i));
        contribution_dark_B_MC{3}(i,i_MC) = decay_dark_B_summer(i)./(decay_pH_B_summer(i) + decay_dark_B_summer(i) + decay_sun_B_summer(i));
        contribution_sun_B_MC{3}(i,i_MC) = decay_sun_B_summer(i)./(decay_pH_B_summer(i) + decay_dark_B_summer(i) + decay_sun_B_summer(i));
    end
    
    for i = 1:n_B_fall
        decay_pH_B_fall(i) = k_pH*10^(data_B_fall(i,2)-14)*teta_pH^(data_B_fall(i,3)-20);
        decay_dark_B_fall(i) = k_dark*teta_dark^(data_B_fall(i,3)-20);
        decay_sun_B_fall(i) = alpha_sun*data_B_fall(i,5)/(data_B_fall(i,4)*d)*(1 - exp(-data_B_fall(i,4)*d));
        
        contribution_pH_B_MC{4}(i,i_MC) = decay_pH_B_fall(i)./(decay_pH_B_fall(i) + decay_dark_B_fall(i) + decay_sun_B_fall(i));
        contribution_dark_B_MC{4}(i,i_MC) = decay_dark_B_fall(i)./(decay_pH_B_fall(i) + decay_dark_B_fall(i) + decay_sun_B_fall(i));
        contribution_sun_B_MC{4}(i,i_MC) = decay_sun_B_fall(i)./(decay_pH_B_fall(i) + decay_dark_B_fall(i) + decay_sun_B_fall(i));
    end
    
    contribution_pH_B_MC{5}(:,i_MC) =... 
        [contribution_pH_B_MC{1}(:,i_MC);...
         contribution_pH_B_MC{2}(:,i_MC);...
         contribution_pH_B_MC{3}(:,i_MC);...
         contribution_pH_B_MC{4}(:,i_MC)];
    contribution_dark_B_MC{5} = ...
        [contribution_dark_B_MC{1}(:,i_MC);...
        contribution_dark_B_MC{2}(:,i_MC);...
        contribution_dark_B_MC{3}(:,i_MC);...
        contribution_dark_B_MC{4}(:,i_MC)];
    contribution_sun_B_MC{5} = ...
        [contribution_sun_B_MC{1}(:,i_MC);...
        contribution_sun_B_MC{2}(:,i_MC);...
        contribution_sun_B_MC{3}(:,i_MC);...
        contribution_sun_B_MC{4}(:,i_MC)];
                     
    i_MC
end
%%

if option_save == 1
    save('./Contribution analysis/yearlyDecayContribution_MonteCarlo.mat',...
        'contribution_pH_A_MC','contribution_dark_A_MC','contribution_sun_A_MC',...
        'contribution_pH_B_MC','contribution_dark_B_MC','contribution_sun_B_MC');
end

