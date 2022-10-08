%% Main_yearly_environmental_data_format

% This script gathers, formats, and draw graphs of the yearly environmental
% data from HRAP A and B to use for the publication Modelling E. coli
% removal during real domestic wastewater treatment in outdoor pilot scale
% High Rate Algal Ponds.

% In this script, the data is first imported, the associated date vectors
% are corrected to suit Matlab format, an NaN data are removed for each
% vector of a given environmental variable.
% Dates vectors are then created from the available data for pH (for which
% temperature is necessarily available too since the sensor was logging
% both variables). 
% Sunlight intensity and light attenuation coefficients were
% computed though linear interpolation on the time vectors constructed.

% Having constructed the data, subset of data are created for each season
% (spring, summer, autumn, and winter) and are formatted for future use in
% boxplots.

% The created data was stored in the Matlab data file
% 'yearly_environmental_data_formatted.mat'.

clear all

%% Import of pH and T data

pH_T_A = xlsread('./Yearly environmental data/HRAP_A_pH_DO_T_yearly data.xlsx');
pH_T_B = xlsread('./Yearly environmental data/HRAP_B_pH_DO_T_yearly data.xlsx');

% Time, pH and temperature are compiled in a single matrix. The time is
% adapted from .xlsx norm to Matlab norm by adding the 693960 coefficient. 
pH_T_A = [pH_T_A(:,1) + 693960 , pH_T_A(:,2) , pH_T_A(:,6)]; 
pH_T_B = [pH_T_B(:,1) + 693960 , pH_T_B(:,2) , pH_T_B(:,6)]; 

% Dates with missing data are removed from the final analysis
[row_A, ~] = find(isnan(pH_T_A));
pH_T_A(row_A,:) = [];

[row_B, ~] = find(isnan(pH_T_B));
pH_T_B(row_B,:) = [];

%% Import of light attenuation

% The same work carried out for pH and temperature is performed for the
% light attenuation coefficient.
sigma = xlsread('./Yearly environmental data/TSS_light attenuation.xlsx');

sigma_A_meas = [sigma(:,1) + 693960 , sigma(:,3)];
sigma_B_meas = [sigma(:,4) + 693960 , sigma(:,6)];

[row_A,~] = find(isnan(sigma_A_meas));
[row_B,~] = find(isnan(sigma_B_meas));

sigma_A_meas(row_A,:) = [];
sigma_B_meas(row_B,:) = [];

%% Sunlight intensity

% The same work carried out for pH and temperature is performed for the
% sunlight intensity.

Hs = xlsread('./Yearly environmental data/NIWA_meteo_data_Palmerston_North.xlsx');
Hs = [Hs(:,2) + 693960,Hs(:,13)];
[row,~] = find(isnan(Hs));
Hs(row,:) = [];

%% Data polishing:
% The objective is here to take all of the data collected for pH and T and
% find the corresponding value for sigma and Hs in both pond. THis value is
% calcualted from the time interpolation of the two data poinds closest
% recorded before and after the data poind of interest.

% First, pH and Temperature are limited to the time interval decided in
% based on the available data for all parameters (see manuscript)

t_A_min = pH_T_A(1,1);
t_A_max = datenum(2017,06,01);

t_B_min = pH_T_B(1,1);
t_B_max = pH_T_B(end,1);

time_A = pH_T_A(:,1);
index_A = t_A_min <= time_A & time_A < t_A_max;
time_A = time_A(index_A);
pH_A = pH_T_A(index_A,2);
T_A =  pH_T_A(index_A,3);

time_B = pH_T_B(:,1);
index_B = t_B_min <= time_B & time_B < t_B_max;
time_B = time_B(index_A);
pH_B = pH_T_B(index_A,2);
T_B =  pH_T_B(index_A,3);


% Hs and sigma are computed on the time vector defined above by linear
% interpolation. To achieve this, the two non-NaN values framing a given
% date point are identified first to perform linear interpolation on these
% values and their corresponding dates.

n_A = length(time_A);
Hs_A = NaN(n_A,1);
sigma_A = NaN(n_A,1);

for i = 1:n_A
    
    k_Hs_A = 1;
    k_sigma_A = 1;
    
    t_A = time_A(i);
    
    if sigma_A_meas(k_sigma_A,1) > t_A
        disp('Error: no prior data for sigma');
    else
        while sigma_A_meas(k_sigma_A,1) <= t_A
            k_sigma_A = k_sigma_A + 1;
        end
        t_sigma_A_0 = sigma_A_meas(k_sigma_A - 1,1);
        sigma_A_0 = sigma_A_meas(k_sigma_A - 1,2);
        t_sigma_A_1 = sigma_A_meas(k_sigma_A,1);
        sigma_A_1 = sigma_A_meas(k_sigma_A,2);
    end
    
    if Hs(k_Hs_A,1) > t_A
        disp('Error: no prior data for Hs');
    else
        while Hs(k_Hs_A,1) <= t_A
            k_Hs_A = k_Hs_A + 1;
        end
        t_Hs_A_0 = Hs(k_Hs_A - 1,1);
        Hs_A_0 = Hs(k_Hs_A - 1,2);
        t_Hs_A_1 = Hs(k_Hs_A,1);
        Hs_A_1 = Hs(k_Hs_A,2);
    end
    
    sigma_A(i) = time_interpolation(t_sigma_A_0,sigma_A_0,t_sigma_A_1,sigma_A_1,t_A);
    Hs_A(i) = time_interpolation(t_Hs_A_0,Hs_A_0,t_Hs_A_1,Hs_A_1,t_A);
end

n_B = length(time_B);
Hs_B = NaN(n_B,1);
sigma_B = NaN(n_B,1);

for i = 1:n_B
    
    k_Hs_B = 1;
    k_sigma_B = 1;
    
    t_A = time_B(i);
    
    if sigma_B_meas(k_sigma_B,1) > t_A
        disp('Error: no prior data for sigma');
    else
        while sigma_B_meas(k_sigma_B,1) <= t_A
            k_sigma_B = k_sigma_B + 1;
        end
        t_sigma_B_0 = sigma_B_meas(k_sigma_B - 1,1);
        sigma_B_0 = sigma_B_meas(k_sigma_B - 1,2);
        t_sigma_B_1 = sigma_B_meas(k_sigma_B,1);
        sigma_B_1 = sigma_B_meas(k_sigma_B,2);
    end
    
    if Hs(k_Hs_B,1) > t_A
        disp('Error: no prior data for Hs');
    else
        while Hs(k_Hs_B,1) <= t_A
            k_Hs_B = k_Hs_B + 1;
        end
        t_Hs_B_0 = Hs(k_Hs_B - 1,1);
        Hs_B_0 = Hs(k_Hs_B - 1,2);
        t_Hs_B_1 = Hs(k_Hs_B,1);
        Hs_B_1 = Hs(k_Hs_B,2);
    end
    
    sigma_B(i) = time_interpolation(t_sigma_B_0,sigma_B_0,t_sigma_B_1,sigma_B_1,t_A);
    Hs_B(i) = time_interpolation(t_Hs_B_0,Hs_B_0,t_Hs_B_1,Hs_B_1,t_A);
end

% The data are compiled in a single matrix per HRAP
data_A = [time_A,pH_A,T_A,sigma_A,Hs_A];
data_B = [time_B,pH_B,T_B,sigma_B,Hs_B];

l_A = size(data_A,1);
l_B = size(data_B,1);


%% Seasonal data separation

% Four data subsets are created to gather data from each season.

t0_winter = datenum(2016,6,1);
t0_spring = datenum(2016,9,1);
t0_summer = datenum(2016,12,1);
t0_fall = datenum(2017,3,1);
t1_fall = datenum(2017,6,1);

data_A_winter = [];
data_A_spring = [];
data_A_summer = [];
data_A_fall = [];

for i = 1:l_A
    if data_A(i,1) >= t0_winter && data_A(i,1) < t0_spring
        data_A_winter = [data_A_winter ; data_A(i,:)];
    end
    if data_A(i,1) >= t0_spring && data_A(i,1) < t0_summer
        data_A_spring = [data_A_spring ; data_A(i,:)];
    end
    if data_A(i,1) >= t0_summer && data_A(i,1) < t0_fall
        data_A_summer = [data_A_summer ; data_A(i,:)];
    end
    if data_A(i,1) >= t0_fall && data_A(i,1) < t1_fall
        data_A_fall = [data_A_fall ; data_A(i,:)];
    end
end


sigma_A_meas_winter = [];
sigma_A_meas_spring = [];
sigma_A_meas_summer = [];
sigma_A_meas_fall = [];

for i = 1:length(sigma_A_meas)
    if sigma_A_meas(i,1) >= t0_winter && sigma_A_meas(i,1) < t0_spring
        sigma_A_meas_winter = [sigma_A_meas_winter ; sigma_A_meas(i,2)];
    end
    if sigma_A_meas(i,1) >= t0_spring && sigma_A_meas(i,1) < t0_summer
        sigma_A_meas_spring = [sigma_A_meas_spring ; sigma_A_meas(i,2)];
    end
    if sigma_A_meas(i,1) >= t0_summer && sigma_A_meas(i,1) < t0_fall
        sigma_A_meas_summer = [sigma_A_meas_summer ; sigma_A_meas(i,2)];
    end
    if sigma_A_meas(i,1) >= t0_fall && sigma_A_meas(i,1) < t1_fall
        sigma_A_meas_fall = [sigma_A_meas_fall ; sigma_A_meas(i,2)];
    end
end

Hs_winter = [];
Hs_spring = [];
Hs_summer = [];
Hs_fall = [];

for i = 1:length(Hs)
    if Hs(i,1) >= t0_winter && Hs(i,1) < t0_spring
        Hs_winter = [Hs_winter ; Hs(i,2)];
    end
    if Hs(i,1) >= t0_spring && Hs(i,1) < t0_summer
        Hs_spring = [Hs_spring ; Hs(i,2)];
    end
    if Hs(i,1) >= t0_summer && Hs(i,1) < t0_fall
        Hs_summer = [Hs_summer ; Hs(i,2)];
    end
    if Hs(i,1) >= t0_fall && Hs(i,1) < t1_fall
        Hs_fall = [Hs_fall ; Hs(i,2)];
    end
end

data_B_winter = [];
data_B_spring = [];
data_B_summer = [];
data_B_fall = [];

for i = 1:l_B
    if data_B(i,1) >= t0_winter && data_B(i,1) < t0_spring
        data_B_winter = [data_B_winter ; data_B(i,:)];
    end
    if data_B(i,1) >= t0_spring && data_B(i,1) < t0_summer
        data_B_spring = [data_B_spring ; data_B(i,:)];
    end
    if data_B(i,1) >= t0_summer && data_B(i,1) < t0_fall
        data_B_summer = [data_B_summer ; data_B(i,:)];
    end
    if data_B(i,1) >= t0_fall && data_B(i,1) < t1_fall
        data_B_fall = [data_B_fall ; data_B(i,:)];
    end
end

sigma_B_meas_winter = [];
sigma_B_meas_spring = [];
sigma_B_meas_summer = [];
sigma_B_meas_fall = [];

for i = 1:length(sigma_B_meas)
    if sigma_B_meas(i,1) >= t0_winter && sigma_B_meas(i,1) < t0_spring
        sigma_B_meas_winter = [sigma_B_meas_winter ; sigma_B_meas(i,2)];
    end
    if sigma_B_meas(i,1) >= t0_spring && sigma_B_meas(i,1) < t0_summer
        sigma_B_meas_spring = [sigma_B_meas_spring ; sigma_B_meas(i,2)];
    end
    if sigma_B_meas(i,1) >= t0_summer && sigma_B_meas(i,1) < t0_fall
        sigma_B_meas_summer = [sigma_B_meas_summer ; sigma_B_meas(i,2)];
    end
    if sigma_B_meas(i,1) >= t0_fall && sigma_B_meas(i,1) < t1_fall
        sigma_B_meas_fall = [sigma_B_meas_fall ; sigma_B_meas(i,2)];
    end
end

l_A_winter = size(data_A_winter,1);
l_A_spring = size(data_A_spring,1);
l_A_summer = size(data_A_summer,1);
l_A_fall = size(data_A_fall,1);

l_B_winter = size(data_B_winter,1);
l_B_spring = size(data_B_spring,1);
l_B_summer = size(data_B_summer,1);
l_B_fall = size(data_B_fall,1);


%% Formatting of the datasets for seasonal boxplots
% Because box plots will be drawn for this study, vectors are created
% to be adapted to drawing such figures

% HRAP A

% For sunlght intensity, 0 values (i.e. night time) are ignored
no_zero = Hs(:,2) ~= 0;
Hs_day = Hs(no_zero,2);
Hs_plot = NaN(length(Hs_day),5);
Hs_plot(:,1) = Hs_day;
Hs_plot(:,2) = [Hs_winter(Hs_winter ~= 0);NaN(length(Hs_day) - length(Hs_winter(Hs_winter ~= 0)),1)];
Hs_plot(:,3) = [Hs_spring(Hs_spring ~= 0);NaN(length(Hs_day) - length(Hs_spring(Hs_spring ~= 0)),1)];
Hs_plot(:,4) = [Hs_summer(Hs_summer ~= 0);NaN(length(Hs_day) - length(Hs_summer(Hs_summer ~= 0)),1)];
Hs_plot(:,5) = [Hs_fall(Hs_fall ~= 0);NaN(length(Hs_day) - length(Hs_fall(Hs_fall ~= 0)),1)];

T_A_plot = NaN(l_A,5);
T_A_plot(:,1) = data_A(:,3);
T_A_plot(:,2) = [data_A_winter(:,3);NaN(l_A - length(data_A_winter(:,3)),1)];
T_A_plot(:,3) = [data_A_spring(:,3);NaN(l_A - length(data_A_spring(:,3)),1)];
T_A_plot(:,4) = [data_A_summer(:,3);NaN(l_A - length(data_A_summer(:,3)),1)];
T_A_plot(:,5) = [data_A_fall(:,3);NaN(l_A - length(data_A_fall(:,3)),1)];

pH_A_plot = NaN(l_A,5);
pH_A_plot(:,1) = data_A(:,2);
pH_A_plot(:,2) = [data_A_winter(:,2);NaN(l_A - length(data_A_winter(:,2)),1)];
pH_A_plot(:,3) = [data_A_spring(:,2);NaN(l_A - length(data_A_spring(:,2)),1)];
pH_A_plot(:,4) = [data_A_summer(:,2);NaN(l_A - length(data_A_summer(:,2)),1)];
pH_A_plot(:,5) = [data_A_fall(:,2);NaN(l_A - length(data_A_fall(:,2)),1)];

sigma_A_plot = NaN(size(sigma_A_meas,1),5);
sigma_A_plot(:,1) = sigma_A_meas(:,2);
sigma_A_plot(:,2) = [sigma_A_meas_winter ; NaN(size(sigma_A_meas,1) - length(sigma_A_meas_winter),1)];
sigma_A_plot(:,3) = [sigma_A_meas_spring ; NaN(size(sigma_A_meas,1) - length(sigma_A_meas_spring),1)];
sigma_A_plot(:,4) = [sigma_A_meas_summer ; NaN(size(sigma_A_meas,1) - length(sigma_A_meas_summer),1)];
sigma_A_plot(:,5) = [sigma_A_meas_fall ; NaN(size(sigma_A_meas,1) - length(sigma_A_meas_fall),1)];


% HRAP B

T_B_plot = NaN(l_B,5);
T_B_plot(:,1) = data_B(:,3);
T_B_plot(:,2) = [data_B_winter(:,3);NaN(l_B - length(data_B_winter(:,3)),1)];
T_B_plot(:,3) = [data_B_spring(:,3);NaN(l_B - length(data_B_spring(:,3)),1)];
T_B_plot(:,4) = [data_B_summer(:,3);NaN(l_B - length(data_B_summer(:,3)),1)];
T_B_plot(:,5) = [data_B_fall(:,3);NaN(l_B - length(data_B_fall(:,3)),1)];

pH_B_plot = NaN(l_B,5);
pH_B_plot(:,1) = data_B(:,2);
pH_B_plot(:,2) = [data_B_winter(:,2);NaN(l_B - length(data_B_winter(:,2)),1)];
pH_B_plot(:,3) = [data_B_spring(:,2);NaN(l_B - length(data_B_spring(:,2)),1)];
pH_B_plot(:,4) = [data_B_summer(:,2);NaN(l_B - length(data_B_summer(:,2)),1)];
pH_B_plot(:,5) = [data_B_fall(:,2);NaN(l_B - length(data_B_fall(:,2)),1)];

sigma_B_plot = NaN(size(sigma_B_meas,1),5);
sigma_B_plot(:,1) = sigma_B_meas(:,2);
sigma_B_plot(:,2) = [sigma_B_meas_winter ; NaN(size(sigma_B_meas,1) - length(sigma_B_meas_winter),1)];
sigma_B_plot(:,3) = [sigma_B_meas_spring ; NaN(size(sigma_B_meas,1) - length(sigma_B_meas_spring),1)];
sigma_B_plot(:,4) = [sigma_B_meas_summer ; NaN(size(sigma_B_meas,1) - length(sigma_B_meas_summer),1)];
sigma_B_plot(:,5) = [sigma_B_meas_fall ; NaN(size(sigma_B_meas,1) - length(sigma_B_meas_fall),1)];



%% Data export
% The data is saved in .mat object named
% 'yearly_environmental_data_formatted.mat'.

save('./Yearly environmental data/yearly_environmental_data_formatted.mat',...
    'data_A'        , 'data_B',...
    'data_A_winter' , 'data_B_winter',...
    'data_A_spring' , 'data_B_spring',...
    'data_A_fall'   , 'data_B_fall',...
    'data_A_summer' , 'data_B_summer',...
    'Hs_plot',...
    'T_A_plot' , 'pH_A_plot' , 'sigma_A_plot',...
    'T_B_plot' , 'pH_B_plot' , 'sigma_B_plot'...
    );





