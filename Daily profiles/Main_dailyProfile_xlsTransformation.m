%% Main_dailyProfile_xlsTransf

% This script aims at transfroming the excel file of raw data from daily
% profile experiment and to transform it into a another excel file that
% will be used for data analysis in relation with daily profiles in the
% publication 'Modelling E. coli removal during real domestic wastewater
% treatment in outdoor pilot scale High Rate Algal Ponds'

% Briefly, time vector are generated on defined time with the desired time
% step. Environmental data from daily profile (Time, pH, temperature and
% sunlight intensity) are computed on the corresponding time stamp by
% linear interpolation of the measured data. E. coli values are only
% replaced on the corresponding time step but are left as NaN on other time
% step.

%% Data import: import of complementary data
% Here, data such as influent flowrate, influent E. coli cell count are
% imported and saved in the final result file.

[~,~,A] = xlsread('/Daily profiles/Daily_profile_rawdata_original.xlsx','Measured parameters');    
xlswrite('/Daily profiles/Daily_profile_transformedData.xlsx',A,'Model','A1')
    
%% Daily profile data import
dp_name = {'30_09' , '12_10' , '29_10' , '17_11' , '03_02' , '10_02', '16_03'};
n_exp = 7;
time_data = cell(n_exp,1);
pH_data = cell(n_exp,1);
DO_data = cell(n_exp,1);
temp_data = cell(n_exp,1);
sun_data = cell(n_exp,1);
coli_data = cell(n_exp,1);
coli_data_IC_minus = cell(n_exp,1);
coli_data_IC_plus = cell(n_exp,1);
n_data = cell(n_exp,1);
new_data = cell(n_exp,1);


for i = 1:n_exp
    
%     i = 1 % just a line to debug that needs to be deleted
    A = xlsread('./Daily profiles/Daily_profile_rawdata_original.xlsx',dp_name{i});
    time_data{i} = A(:,1);
    sun_data{i} = A(:,5);
    pH_data{i} = A(:,2);
    DO_data{i} = A(:,4);
    temp_data{i} = A(:,3);
    coli_data{i} = A(:,6);
    coli_data_IC_minus{i} = A(:,7);
    coli_data_IC_plus{i} = A(:,8);
    
    
    n_data{i} = length(time_data{i});
    
    t_start = time_data{i}(1) + 693960;
    t_end = time_data{i}(end) + 693960;
    t_start = floor(t_start*24*60)/24/60;
    t_end = floor(t_end*24*60)/24/60;
    
    % Initialisation: a number of point equivalent to a 1 s time step for
    % the time data is defined and corresponding NaN vectors for outputs
    % are created
    new_time_data = ((t_start*24*60):(t_end*24*60))/24/60;
    new_n = length(new_time_data);
    new_sun_data = NaN(new_n,1);
    new_pH_data = NaN(new_n,1);
    new_DO_data = NaN(new_n,1);
    new_temp_data = NaN(new_n,1);
    new_coli_data = NaN(new_n,1);
    new_coli_data_IC_minus = NaN(new_n,1);
    new_coli_data_IC_plus = NaN(new_n,1);
    
    % Initialisation of new vectors
    new_sun_data(1) = sun_data{i}(1);
    new_pH_data(1) = pH_data{i}(1);
    new_temp_data(1) = temp_data{i}(1);
    new_DO_data(1) = DO_data{i}(1);
    new_coli_data(1) = coli_data{i}(1);
    new_coli_data_IC_minus(1) = coli_data_IC_minus{i}(1);
    new_coli_data_IC_plus(1) = coli_data_IC_plus{i}(1);
    
    old_p = 2;
    new_p = 2;
    while new_p <= new_n && old_p <= n_data{i}
        local_p = new_p-1;
        while new_time_data(new_p) < time_data{i}(old_p)+693960 && new_p < new_n
            new_p = new_p + 1;
        end
        new_sun_data(new_p) = sun_data{i}(old_p);
        new_pH_data(new_p) = pH_data{i}(old_p);
        new_DO_data(new_p) = DO_data{i}(old_p);
        new_temp_data(new_p) = temp_data{i}(old_p);
        
        % Lieaer interpolation on the data between local_p + 1 and new_p -
        % 1
        for new_k = local_p + 1:new_p-1
            new_sun_data(new_k) = new_sun_data(local_p)+(new_sun_data(new_p) - new_sun_data(local_p))/(new_time_data(new_p)-new_time_data(local_p))*(new_time_data(new_k)-new_time_data(local_p));
            new_pH_data(new_k) = new_pH_data(local_p)+(new_pH_data(new_p) - new_pH_data(local_p))/(new_time_data(new_p)-new_time_data(local_p))*(new_time_data(new_k)-new_time_data(local_p));
            new_temp_data(new_k) = new_temp_data(local_p)+(new_temp_data(new_p) - new_temp_data(local_p))/(new_time_data(new_p)-new_time_data(local_p))*(new_time_data(new_k)-new_time_data(local_p));
            new_DO_data(new_k) = new_DO_data(local_p)+(new_DO_data(new_p) - new_DO_data(local_p))/(new_time_data(new_p)-new_time_data(local_p))*(new_time_data(new_k)-new_time_data(local_p));
        end
    
        % E. coli cell count
        if isnan(coli_data{i}(old_p)) == 0
            new_coli_data(new_p) = coli_data{i}(old_p);
            new_coli_data_IC_minus(new_p) = coli_data_IC_minus{i}(old_p);
            new_coli_data_IC_plus(new_p) = coli_data_IC_plus{i}(old_p);
        end
        
        new_p = new_p + 1;
        old_p = old_p + 1;
    end
    
    new_coli_data(end) = coli_data{i}(end);
    
    B = NaN(new_n,6);
    B(:,1) = new_time_data;
    B(:,2) = new_pH_data;
    B(:,3) = new_temp_data;
    B(:,4) = new_DO_data;
    B(:,5) = new_sun_data;
    B(:,6) = new_coli_data;
    B(:,7) = new_coli_data_IC_minus;
    B(:,8) = new_coli_data_IC_plus;
    
    % Add xlswrite when ready
    xlswrite('./Daily profiles/Daily_profile_transformedData.xlsx',B,dp_name{i},'A2')
    
end
    