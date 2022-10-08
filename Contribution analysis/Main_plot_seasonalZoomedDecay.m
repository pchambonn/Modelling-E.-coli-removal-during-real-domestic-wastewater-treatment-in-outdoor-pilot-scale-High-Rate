%% Main_plot_seasonalZoomedDecay

% This script enables to plot Figure 4 and Figure S10.1 of the publication
% 'Modelling E. coli removal during real domestic wastewater treatment in
% outdoor pilot scale High Rate Algal Ponds'.

%%
clear all

load('./Contribution analysis/yearlyDecayValue.mat')
load('./Yearly environmental data/yearly_environmental_data_formatted.mat')

%% Seasonally zoomed in plots

% HRAP A
figure(12), clf, hold on
fs = 16; 


subplot(2,2,1)

decay_timeplot_winter = [decay_dark_A{1},decay_pH_A{1},decay_sun_A{1}];
time_data_loc = datetime(0,1,floor(data_A_winter(:,1)),0,0,3600*24*(data_A_winter(:,1) -floor(data_A_winter(:,1))));
    
 
b_winter = bar(time_data_loc,decay_timeplot_winter,85,'stacked');

b_winter(1).FaceColor = [0.2422 0.1504 0.6603];
b_winter(2).FaceColor = [0.0704 0.7457 0.7258];
b_winter(3).FaceColor = [0.9769 0.9839 0.0805];


ax1 = gca;
ax1.YLim = [0 50];
x0 = find(time_data_loc >= datetime(2016,07,23,0,0,0) & time_data_loc <= datetime(2016,07,23,0,15,0)); 
x1 = find(time_data_loc >= datetime(2016,07,30,0,0,0) & time_data_loc <= datetime(2016,07,30,0,15,0)); 
ax1.FontSize = fs - 4;
legend(b_winter,{'Uncharacterized dark decay','pH toxicity','Sunlight direct damage'},'FontSize',fs - 6,'Location','North')
datetick('x','dd-mmm-yyyy')
ax1.XTick = [time_data_loc(x0,1) time_data_loc(x1,1)];
ax1.XTickLabel = {datestr(time_data_loc(x0,1),1) , datestr(time_data_loc(x1,1),1)};
ax1.XLim = [time_data_loc(x0,1) time_data_loc(x1,1)];

title('Winter','FontSize',fs)
ylabel('Mechanism decay rate ({\cdot}d^-^1)','FontSize',fs-3,'FontWeight','bold')


subplot(2,2,2)

decay_timeplot_spring = ...
    [decay_dark_A{2},decay_pH_A{2},decay_sun_A{2}];
time_data_loc = datetime(0,1,floor(data_A_spring(:,1)),0,0,3600*24*(data_A_spring(:,1) -floor(data_A_spring(:,1))));
    
b_spring = bar(time_data_loc,decay_timeplot_spring,900,'stacked');

b_spring(1).FaceColor = [0.2422 0.1504 0.6603];
b_spring(2).FaceColor = [0.0704 0.7457 0.7258];
b_spring(3).FaceColor = [0.9769 0.9839 0.0805];

ax1 = gca;
ax1.YLim = [0 50];
x0 = min(find(time_data_loc >= datetime(2016,10,02,0,0,0) & time_data_loc <= datetime(2016,10,02,0,15,0))); 
x1 = min(find(time_data_loc >= datetime(2016,10,09,0,0,0) & time_data_loc <= datetime(2016,10,09,0,15,0))); 
ax1.FontSize = fs - 4;
datetick('x','dd-mmm-yyyy')
ax1.XTick = [time_data_loc(x0,1) time_data_loc(x1,1)];
ax1.XTickLabel = {datestr(time_data_loc(x0,1),1) , datestr(time_data_loc(x1,1),1)};
ax1.XLim = [time_data_loc(x0,1) time_data_loc(x1,1)];

title('Spring','FontSize',fs)


subplot(2,2,4)

decay_timeplot_summer = [decay_dark_A{3},decay_pH_A{3},decay_sun_A{3}];
time_data_loc = datetime(0,1,floor(data_A_summer(:,1)),0,0,3600*24*(data_A_summer(:,1) -floor(data_A_summer(:,1))));
     
b_summer = bar(time_data_loc,decay_timeplot_summer,4,'stacked');

b_summer(1).FaceColor = [0.2422 0.1504 0.6603];
b_summer(2).FaceColor = [0.0704 0.7457 0.7258];
b_summer(3).FaceColor = [0.9769 0.9839 0.0805];


ax1 = gca;
ax1.YLim = [0 85];
x0 = min(find(time_data_loc >= datetime(2017,02,12,0,0,0) & time_data_loc <= datetime(2017,02,12,0,15,0))); 
x1 = min(find(time_data_loc >= datetime(2017,02,19,0,0,0) & time_data_loc <= datetime(2017,02,19,0,15,0))); 
ax1.FontSize = fs - 4;
datetick('x','dd-mmm-yyyy')
ax1.XTick = [time_data_loc(x0,1) time_data_loc(x1,1)];
ax1.XTickLabel = {datestr(time_data_loc(x0,1),1) , datestr(time_data_loc(x1,1),1)};
ax1.XLim = [time_data_loc(x0,1) time_data_loc(x1,1)];

title('Summer','FontSize',fs)


subplot(2,2,3)

decay_timeplot_fall = [decay_dark_A{4},decay_pH_A{4},decay_sun_A{4}];
time_data_loc = datetime(0,1,floor(data_A_fall(:,1)),0,0,3600*24*(data_A_fall(:,1) -floor(data_A_fall(:,1))));
     
b_fall = bar(time_data_loc,decay_timeplot_fall,10,'stacked');

b_fall(1).FaceColor = [0.2422 0.1504 0.6603];
b_fall(2).FaceColor = [0.0704 0.7457 0.7258];
b_fall(3).FaceColor = [0.9769 0.9839 0.0805];

ax1 = gca;
ax1.YLim = [0 50];
x0 = min(find(time_data_loc >= datetime(2017,04,09,0,0,0) & time_data_loc <= datetime(2017,04,09,0,15,0))); 
x1 = min(find(time_data_loc >= datetime(2017,04,16,0,0,0) & time_data_loc <= datetime(2017,04,16,0,15,0))); 
ax1.FontSize = fs - 4;
datetick('x','dd-mmm-yyyy')
ax1.XTick = [time_data_loc(x0,1) time_data_loc(x1,1)];
ax1.XTickLabel = {datestr(time_data_loc(x0,1),1) , datestr(time_data_loc(x1,1),1)};
ax1.XLim = [time_data_loc(x0,1) time_data_loc(x1,1)];

ylabel('Mechanism decay rate ({\cdot}d^-^1)','FontSize',fs-3,'FontWeight','bold')
title('Fall','FontSize',fs)

fig = gcf;
fig.Position = [40 40 0.9*900 0.9*800];




%%
% HRAP B
figure(22), clf, hold on
fs = 16; 

 
subplot(2,2,1)

decay_timeplot_winter = [decay_dark_B{1},decay_pH_B{1},decay_sun_B{1}];
time_data_loc = datetime(0,1,floor(data_B_winter(:,1)),0,0,3600*24*(data_B_winter(:,1) -floor(data_B_winter(:,1))));
    
 
b_winter = bar(time_data_loc,decay_timeplot_winter,2,'stacked');

b_winter(1).FaceColor = [0.2422 0.1504 0.6603];
b_winter(2).FaceColor = [0.0704 0.7457 0.7258];
b_winter(3).FaceColor = [0.9769 0.9839 0.0805];

ax1 = gca;
ax1.FontSize = fs - 4;
ax1.YLim = [0 50];
legend(b_winter,{'Uncharacterized dark decay','pH toxicity','Sunlight direct damage'},'FontSize',fs - 6,'Location','north')
x0 = min(find(time_data_loc >= datetime(2016,08,17,0,0,0) & time_data_loc <= datetime(2016,08,17,0,15,0))); 
x1 = min(find(time_data_loc >= datetime(2016,08,24,0,0,0) & time_data_loc <= datetime(2016,08,24,0,15,0))); 
ax1.FontSize = fs - 4;
datetick('x','dd-mmm-yyyy')
ax1.XTick = [time_data_loc(x0,1) time_data_loc(x1,1)];
ax1.XTickLabel = {datestr(time_data_loc(x0,1),1) , datestr(time_data_loc(x1,1),1)};
ax1.XLim = [time_data_loc(x0,1) time_data_loc(x1,1)];

title('Winter','FontSize',fs)
ylabel('Mechanism decay rate ({\cdot}d^-^1)','FontSize',fs-3,'FontWeight','bold')


subplot(2,2,2)

decay_timeplot_spring = ...
    [decay_dark_B{2},decay_pH_B{2},decay_sun_B{2}];
time_data_loc = datetime(0,1,floor(data_B_spring(:,1)),0,0,3600*24*(data_B_spring(:,1) -floor(data_B_spring(:,1))));
    
 b_spring = bar(time_data_loc,decay_timeplot_spring,100,'stacked');

b_spring(1).FaceColor = [0.2422 0.1504 0.6603];
b_spring(2).FaceColor = [0.0704 0.7457 0.7258];
b_spring(3).FaceColor = [0.9769 0.9839 0.0805];

ax1 = gca;
ax1.YLim = [0 50];
ax1.FontSize = fs - 4;
x0 = min(find(time_data_loc >= datetime(2016,11,23,0,0,0) & time_data_loc <= datetime(2016,11,23,0,15,0))); 
x1 = min(find(time_data_loc >= datetime(2016,11,30,0,0,0) & time_data_loc <= datetime(2016,11,30,0,15,0))); 
ax1.FontSize = fs - 4;
datetick('x','dd-mmm-yyyy')
ax1.XTick = [time_data_loc(x0,1) time_data_loc(x1,1)];
ax1.XTickLabel = {datestr(time_data_loc(x0,1),1) , datestr(time_data_loc(x1,1),1)};
ax1.XLim = [time_data_loc(x0,1) time_data_loc(x1,1)];

title('Spring','FontSize',fs)


subplot(2,2,4)

decay_timeplot_summer = ...
    [decay_dark_B{3},decay_pH_B{3},decay_sun_B{3}];
time_data_loc = datetime(0,1,floor(data_B_summer(:,1)),0,0,3600*24*(data_B_summer(:,1) -floor(data_B_summer(:,1))));
     
b_summer = bar(time_data_loc,decay_timeplot_summer,30,'stacked');

b_summer(1).FaceColor = [0.2422 0.1504 0.6603];
b_summer(2).FaceColor = [0.0704 0.7457 0.7258];
b_summer(3).FaceColor = [0.9769 0.9839 0.0805];

ax1 = gca;
% ax1.YLim = [0 50];
ax1.FontSize = fs - 4;
x0 = min(find(time_data_loc >= datetime(2017,02,11,0,0,0) & time_data_loc <= datetime(2017,02,11,0,15,0))); 
x1 = min(find(time_data_loc >= datetime(2017,02,18,0,0,0) & time_data_loc <= datetime(2017,02,18,0,15,0))); 
ax1.FontSize = fs - 4;
datetick('x','dd-mmm-yyyy')
ax1.XTick = [time_data_loc(x0,1) time_data_loc(x1,1)];
ax1.XTickLabel = {datestr(time_data_loc(x0,1),1) , datestr(time_data_loc(x1,1),1)};
ax1.XLim = [time_data_loc(x0,1) time_data_loc(x1,1)];

title('Summer','FontSize',fs)


subplot(2,2,3)

decay_timeplot_fall = [decay_dark_B{4},decay_pH_B{4},decay_sun_B{4}];
time_data_loc = datetime(0,1,floor(data_B_fall(:,1)),0,0,3600*24*(data_B_fall(:,1) -floor(data_B_fall(:,1))));
     
b_fall = bar(time_data_loc,decay_timeplot_fall,10,'stacked');

b_fall(1).FaceColor = [0.2422 0.1504 0.6603];
b_fall(2).FaceColor = [0.0704 0.7457 0.7258];
b_fall(3).FaceColor = [0.9769 0.9839 0.0805];

ax1 = gca;
ax1.YLim = [0 50];
ax1.FontSize = fs - 4;
x0 = min(find(time_data_loc >= datetime(2017,04,02,0,0,0) & time_data_loc <= datetime(2017,04,02,0,15,0))); 
x1 = min(find(time_data_loc >= datetime(2017,04,09,0,0,0) & time_data_loc <= datetime(2017,04,09,0,15,0))); 
ax1.FontSize = fs - 4;
datetick('x','dd-mmm-yyyy')
ax1.XTick = [time_data_loc(x0,1) time_data_loc(x1,1)];
ax1.XTickLabel = {datestr(time_data_loc(x0,1),1) , datestr(time_data_loc(x1,1),1)};
ax1.XLim = [time_data_loc(x0,1) time_data_loc(x1,1)];

title('Fall','FontSize',fs)
ylabel('Mechanism decay rate ({\cdot}d^-^1)','FontSize',fs-3,'FontWeight','bold')


fig = gcf;
fig.Position = [40 40 0.9*900 0.9*800];

