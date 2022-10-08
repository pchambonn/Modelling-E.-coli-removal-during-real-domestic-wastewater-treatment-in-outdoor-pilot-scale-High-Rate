%% Main_figure_daily_profile_measurements

% This script aims at drawing the figure describing measurements from
% daily profiles of sunlight intensity, pH, temperature, and E. coli cell
% counts. This figure corresponds to Figure S1.1 in the publication
% PREDICTION OF E. COLI DECAY MECHANISMS DURING WASTEWATER TREATMENT IN
% HIGH RATE ALGAL PONDS.

clear all

%% pH, temperature, sunlight
% The environmental and E. coli data for daily profile is here imported
% from the relevant Excel data file. A data matrix for each day of the
% seven days available is stored in one dimension of a seven dimention cell
% array.
raw_data = cell(7,1);

for i = 2:8
    raw_data{i-1} = xlsread('/Daily profiles/Daily_profile_rawdata_original.xlsx',i);
end

%%
% Figure S1.1 is now drawn pH, temperature, sunlight intensity, and E. coli
% cell count occupying each one quadrant.  using the same symbol code for a
% given day in each quadrant.

figure(1), clf, hold on

fs = 13;

subplot(2,2,1), hold on

    plot(raw_data{1}(:,1)-floor(raw_data{1}(1,1)),raw_data{1}(:,2),'kx')
    plot(raw_data{2}(:,1)-floor(raw_data{2}(1,1)),raw_data{2}(:,2),'ko')
    plot(raw_data{3}(:,1)-floor(raw_data{3}(1,1)),raw_data{3}(:,2),'k*')
    plot(raw_data{4}(:,1)-floor(raw_data{4}(1,1)),raw_data{4}(:,2),'k.')
    plot(raw_data{5}(:,1)-floor(raw_data{5}(1,1)),raw_data{5}(:,2),'k+')
    plot(raw_data{6}(:,1)-floor(raw_data{6}(1,1)),raw_data{6}(:,2),'ks')
    plot(raw_data{7}(:,1)-floor(raw_data{7}(1,1)),raw_data{7}(:,2),'kp')
    
ax = gca;
ax.XLim = [9/24 1 + 9/24];
datetick('x','HH:MM');
ax.FontSize = fs - 2;
ylabel('pH','FontSize',fs,'FontWeight','bold')
ax.XTickLabelRotation = 45

subplot(2,2,2), hold on

    plot(raw_data{1}(:,1)-floor(raw_data{1}(1,1)),raw_data{1}(:,3),'kx')
    plot(raw_data{2}(:,1)-floor(raw_data{2}(1,1)),raw_data{2}(:,3),'ko')
    plot(raw_data{3}(:,1)-floor(raw_data{3}(1,1)),raw_data{3}(:,3),'k*')
    plot(raw_data{4}(:,1)-floor(raw_data{4}(1,1)),raw_data{4}(:,3),'k.')
    plot(raw_data{5}(:,1)-floor(raw_data{5}(1,1)),raw_data{5}(:,3),'k+')
    plot(raw_data{6}(:,1)-floor(raw_data{6}(1,1)),raw_data{6}(:,3),'ks')
    plot(raw_data{7}(:,1)-floor(raw_data{7}(1,1)),raw_data{7}(:,3),'kp')
    
ax = gca;
ax.XLim = [9/24 1 + 9/24]; ;
datetick('x','HH:MM');
ax.FontSize = fs - 2;
ylabel('Temperature (°C)','FontSize',fs,'FontWeight','bold')
ax.XTickLabelRotation = 45

subplot(2,2,3), hold on

    plot(raw_data{1}(1:4:end,1)-floor(raw_data{1}(1,1)),raw_data{1}(1:4:end,5),'kx--')
    plot(raw_data{2}(1:4:end,1)-floor(raw_data{2}(1,1)),raw_data{2}(1:4:end,5),'ko--')
    plot(raw_data{3}(1:4:end,1)-floor(raw_data{3}(1,1)),raw_data{3}(1:4:end,5),'k*--')
    plot(raw_data{4}(1:4:end,1)-floor(raw_data{4}(1,1)),raw_data{4}(1:4:end,5),'k.--')
    plot(raw_data{5}(1:4:end,1)-floor(raw_data{5}(1,1)),raw_data{5}(1:4:end,5),'k+--')
    plot(raw_data{6}(1:4:end,1)-floor(raw_data{6}(1,1)),raw_data{6}(1:4:end,5),'ks--')
    plot(raw_data{7}(1:4:end,1)-floor(raw_data{7}(1,1)),raw_data{7}(1:4:end,5),'kp--')

ax = gca;
ax.XLim = [9/24 1 + 9/24]; ;
datetick('x','HH:MM');
ax.FontSize = fs - 2;
ylabel('Sunlight intensity (W.m^{-2})','FontSize',fs,'FontWeight','bold')
xlabel('Time (HH:MM)','FontSize',fs,'FontWeight','bold')
ax.XTickLabelRotation = 45
legend('30 Sep - 01 Oct 2015','12 - 13 Oct 2015', '28 - 29 0ct 2015','16 - 17 Nov 2015','03 - 04 Feb 2016' ,'10 - 11 Feb 2016', '16 - 17 Mar 2016','FontSize',fs-5)

subplot(2,2,4), hold on
coli = cell(7,1);
for i = 1:7
    coli{i} = fillmissing(raw_data{i}(:,6),'linear');
end
    
    plot(raw_data{1}(:,1)-floor(raw_data{1}(1,1)),log10(raw_data{1}(:,6)),'kx')
    plot(raw_data{1}(:,1)-floor(raw_data{1}(1,1)),log10(coli{1}),'k--')
    plot(raw_data{2}(:,1)-floor(raw_data{2}(1,1)),log10(raw_data{2}(:,6)),'ko')
    plot(raw_data{2}(:,1)-floor(raw_data{2}(1,1)),log10(coli{2}),'k--')
    plot(raw_data{3}(:,1)-floor(raw_data{3}(1,1)),log10(raw_data{3}(:,6)),'k*')
    plot(raw_data{3}(:,1)-floor(raw_data{3}(1,1)),log10(coli{3}),'k--')
    plot(raw_data{4}(:,1)-floor(raw_data{4}(1,1)),log10(raw_data{4}(:,6)),'k.')
    plot(raw_data{4}(:,1)-floor(raw_data{4}(1,1)),log10(coli{4}),'k--')
    plot(raw_data{5}(:,1)-floor(raw_data{5}(1,1)),log10(raw_data{5}(:,6)),'k+')
    plot(raw_data{5}(:,1)-floor(raw_data{5}(1,1)),log10(coli{5}),'k--')
    plot(raw_data{6}(:,1)-floor(raw_data{6}(1,1)),log10(raw_data{6}(:,6)),'ks')
    plot(raw_data{6}(:,1)-floor(raw_data{6}(1,1)),log10(coli{6}),'k--')
    plot(raw_data{7}(:,1)-floor(raw_data{7}(1,1)),log10(raw_data{7}(:,6)),'kp')
    plot(raw_data{7}(:,1)-floor(raw_data{7}(1,1)),log10(coli{7}),'k--')
    

ax = gca;
ax.XLim = [9/24 1 + 9/24];
yticks(3:1:6)
datetick('x','HH:MM');
ax.FontSize = fs - 2;
ylabel({'{\itE. coli} cell count' ; '(log_{10} MPN.100 mL^{-1})'},'FontSize',fs,'FontWeight','bold')
xlabel('Time (HH:MM)','FontSize',fs,'FontWeight','bold')
ax.XTickLabelRotation = 45

fig = gcf;
fig.Position = [50 50 1.05*750 1.05*700];