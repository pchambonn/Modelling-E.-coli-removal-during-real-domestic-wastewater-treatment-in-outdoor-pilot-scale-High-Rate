%% Main_plot_yearlySeasonalContribution

% This script enables to plot Figure 3 and Figure S9.2 of the publication
% 'Modelling E. coli removal during real domestic wastewater treatment in
% outdoor pilot scale High Rate Algal Ponds'.

%%
clear all

load('./Contribution analysis/yearlyDecayContribution.mat')
load('./Daily profiles/inputData_fit_MonteCarlo.mat')
load('./Contribution analysis/yearlyDecayContribution_MonteCarlo.mat')

N_MC = size(contribution_dark_A_MC{1,1},2);

%% Figure S9.2

fs = 18;

figure(1), clf, hold on
% sgtitle('HRAP B','FontSize',fs+2,'FontWeight','bold')
X = [0.75  2.5  4.25];
h1 = boxplot([contribution_dark_A{5} , contribution_pH_A{5} , contribution_sun_A{5}] *100,...
    'positions',X - 0.25,'colors','k','width',0.45,'PlotStyle','traditional','MedianStyle',...
    'line','OutlierSize',3,'whisker',0.7193);
h2 = boxplot([contribution_dark_B{5} , contribution_pH_B{5} , contribution_sun_B{5}] *100,...
    'positions',X + 0.25,'colors',[0.5 0.5 0.5],'width',0.45,'PlotStyle','traditional','MedianStyle',...
    'line','OutlierSize',3,'whisker',0.7193);
ylabel('Mechanism relative contribution (%)','FontSize',fs)
set(gca,'xtick',X,'xticklabel',{'Dark decay' 'pH toxicity' 'Sunlight direct damage'},'FontSize',fs-2,'FontWeight','bold')
ylim([0 100])
xlim([0 5])
set(gca,'LineWidth',3)
set(h1,'LineWidth',2)
set(h2,'LineWidth',2)
yticks(0:20:100)

%Legend work around
h = findobj(gca,'Tag','Box');
hLegend = legend(h([4,1]), {'HRAP A','HRAP B'},'FontSize',fs,'Location','North');

fig = gcf;
fig.Position = [40,40,800,600];


%% Plot from Monte Carlo

% Base value of contribution

total_contribution_dark_A_winter = mean(contribution_dark_A{1});
total_contribution_dark_A_spring = mean(contribution_dark_A{2});
total_contribution_dark_A_summer = mean(contribution_dark_A{3});
total_contribution_dark_A_fall = mean(contribution_dark_A{4});
total_contribution_dark_A_year = mean([total_contribution_dark_A_winter,...
    total_contribution_dark_A_spring,total_contribution_dark_A_summer,...
    total_contribution_dark_A_fall]);

total_contribution_pH_A_winter = mean(contribution_pH_A{1});
total_contribution_pH_A_spring = mean(contribution_pH_A{2});
total_contribution_pH_A_summer = mean(contribution_pH_A{3});
total_contribution_pH_A_fall = mean(contribution_pH_A{4});
total_contribution_pH_A_year = mean([total_contribution_pH_A_winter,...
    total_contribution_pH_A_spring,total_contribution_pH_A_summer,...
    total_contribution_pH_A_fall]);

total_contribution_sun_A_winter = mean(contribution_sun_A{1});
total_contribution_sun_A_spring = mean(contribution_sun_A{2});
total_contribution_sun_A_summer = mean(contribution_sun_A{3});
total_contribution_sun_A_fall = mean(contribution_sun_A{4});
total_contribution_sun_A_year = mean([total_contribution_sun_A_winter,...
    total_contribution_sun_A_spring,total_contribution_sun_A_summer,...
    total_contribution_sun_A_fall]);


% Monte carlo table of mean contribution

total_contribution_dark_A_winter_MC = NaN(N_MC,1);
total_contribution_dark_A_spring_MC = NaN(N_MC,1);
total_contribution_dark_A_summer_MC = NaN(N_MC,1);
total_contribution_dark_A_fall_MC = NaN(N_MC,1);
total_contribution_dark_A_year_MC = NaN(N_MC,1);

total_contribution_pH_A_winter_MC = NaN(N_MC,1);
total_contribution_pH_A_spring_MC = NaN(N_MC,1);
total_contribution_pH_A_summer_MC = NaN(N_MC,1);
total_contribution_pH_A_fall_MC = NaN(N_MC,1);
total_contribution_pH_A_year_MC = NaN(N_MC,1);

total_contribution_sun_A_winter_MC = NaN(N_MC,1);
total_contribution_sun_A_spring_MC = NaN(N_MC,1);
total_contribution_sun_A_summer_MC = NaN(N_MC,1);
total_contribution_sun_A_fall_MC = NaN(N_MC,1);
total_contribution_sun_A_year_MC = NaN(N_MC,1);

for i_MC = 1:N_MC

    total_contribution_dark_A_winter_MC(i_MC) = mean(contribution_dark_A_MC{1}(:,i_MC));
    total_contribution_dark_A_spring_MC(i_MC) = mean(contribution_dark_A_MC{2}(:,i_MC));
    total_contribution_dark_A_summer_MC(i_MC) = mean(contribution_dark_A_MC{3}(:,i_MC));
    total_contribution_dark_A_fall_MC(i_MC) = mean(contribution_dark_A_MC{4}(:,i_MC));
    
    total_contribution_pH_A_winter_MC(i_MC) = mean(contribution_pH_A_MC{1}(:,i_MC));
    total_contribution_pH_A_spring_MC(i_MC) = mean(contribution_pH_A_MC{2}(:,i_MC));
    total_contribution_pH_A_summer_MC(i_MC) = mean(contribution_pH_A_MC{3}(:,i_MC));
    total_contribution_pH_A_fall_MC(i_MC) = mean(contribution_pH_A_MC{4}(:,i_MC));
    
    total_contribution_sun_A_winter_MC(i_MC) = mean(contribution_sun_A_MC{1}(:,i_MC));
    total_contribution_sun_A_spring_MC(i_MC) = mean(contribution_sun_A_MC{2}(:,i_MC));
    total_contribution_sun_A_summer_MC(i_MC) = mean(contribution_sun_A_MC{3}(:,i_MC));
    total_contribution_sun_A_fall_MC(i_MC) = mean(contribution_sun_A_MC{4}(:,i_MC));
end

total_contribution_dark_A_year_MC = (total_contribution_dark_A_winter_MC + total_contribution_dark_A_spring_MC +...
                                     total_contribution_dark_A_summer_MC + total_contribution_dark_A_fall_MC)./4;
                                 
total_contribution_pH_A_year_MC = (total_contribution_pH_A_winter_MC + total_contribution_pH_A_spring_MC +...
                                     total_contribution_pH_A_summer_MC + total_contribution_pH_A_fall_MC)./4;
                                 
total_contribution_sun_A_year_MC = (total_contribution_sun_A_winter_MC + total_contribution_sun_A_spring_MC +...
                                     total_contribution_sun_A_summer_MC + total_contribution_sun_A_fall_MC)./4;

contribution_all_A_winter_MC = [total_contribution_dark_A_winter_MC , total_contribution_pH_A_winter_MC , total_contribution_sun_A_winter_MC];
contribution_all_A_spring_MC = [total_contribution_dark_A_spring_MC , total_contribution_pH_A_spring_MC , total_contribution_sun_A_spring_MC];
contribution_all_A_summer_MC = [total_contribution_dark_A_summer_MC , total_contribution_pH_A_summer_MC , total_contribution_sun_A_summer_MC];
contribution_all_A_fall_MC = [total_contribution_dark_A_fall_MC , total_contribution_pH_A_fall_MC , total_contribution_sun_A_fall_MC];

% HRAP B
total_contribution_dark_B_winter_MC = NaN(N_MC,1);
total_contribution_dark_B_spring_MC = NaN(N_MC,1);
total_contribution_dark_B_summer_MC = NaN(N_MC,1);
total_contribution_dark_B_fall_MC = NaN(N_MC,1);
total_contribution_dark_B_year_MC = NaN(N_MC,1);

total_contribution_pH_B_winter_MC = NaN(N_MC,1);
total_contribution_pH_B_spring_MC = NaN(N_MC,1);
total_contribution_pH_B_summer_MC = NaN(N_MC,1);
total_contribution_pH_B_fall_MC = NaN(N_MC,1);
total_contribution_pH_B_year_MC = NaN(N_MC,1);

total_contribution_sun_B_winter_MC = NaN(N_MC,1);
total_contribution_sun_B_spring_MC = NaN(N_MC,1);
total_contribution_sun_B_summer_MC = NaN(N_MC,1);
total_contribution_sun_B_fall_MC = NaN(N_MC,1);
total_contribution_sun_B_year_MC = NaN(N_MC,1);

for i_MC = 1:N_MC

    total_contribution_dark_B_winter_MC(i_MC) = mean(contribution_dark_B_MC{1}(:,i_MC));
    total_contribution_dark_B_spring_MC(i_MC) = mean(contribution_dark_B_MC{2}(:,i_MC));
    total_contribution_dark_B_summer_MC(i_MC) = mean(contribution_dark_B_MC{3}(:,i_MC));
    total_contribution_dark_B_fall_MC(i_MC) = mean(contribution_dark_B_MC{4}(:,i_MC));

    total_contribution_pH_B_winter_MC(i_MC) = mean(contribution_pH_B_MC{1}(:,i_MC));
    total_contribution_pH_B_spring_MC(i_MC) = mean(contribution_pH_B_MC{2}(:,i_MC));
    total_contribution_pH_B_summer_MC(i_MC) = mean(contribution_pH_B_MC{3}(:,i_MC));
    total_contribution_pH_B_fall_MC(i_MC) = mean(contribution_pH_B_MC{4}(:,i_MC));

    total_contribution_sun_B_winter_MC(i_MC) = mean(contribution_sun_B_MC{1}(:,i_MC));
    total_contribution_sun_B_spring_MC(i_MC) = mean(contribution_sun_B_MC{2}(:,i_MC));
    total_contribution_sun_B_summer_MC(i_MC) = mean(contribution_sun_B_MC{3}(:,i_MC));
    total_contribution_sun_B_fall_MC(i_MC) = mean(contribution_sun_B_MC{4}(:,i_MC));
end

total_contribution_dark_B_year_MC = (total_contribution_dark_B_winter_MC + total_contribution_dark_B_spring_MC +...
                                     total_contribution_dark_B_summer_MC + total_contribution_dark_B_fall_MC)./4;
                                 
total_contribution_pH_B_year_MC = (total_contribution_pH_B_winter_MC + total_contribution_pH_B_spring_MC +...
                                     total_contribution_pH_B_summer_MC + total_contribution_pH_B_fall_MC)./4;
                                 
total_contribution_sun_B_year_MC = (total_contribution_sun_B_winter_MC + total_contribution_sun_B_spring_MC +...
                                     total_contribution_sun_B_summer_MC + total_contribution_sun_B_fall_MC)./4;
                                 

contribution_all_B_winter_MC = [total_contribution_dark_B_winter_MC , total_contribution_pH_B_winter_MC , total_contribution_sun_B_winter_MC];
contribution_all_B_spring_MC = [total_contribution_dark_B_spring_MC , total_contribution_pH_B_spring_MC , total_contribution_sun_B_spring_MC];
contribution_all_B_summer_MC = [total_contribution_dark_B_summer_MC , total_contribution_pH_B_summer_MC , total_contribution_sun_B_summer_MC];
contribution_all_B_fall_MC = [total_contribution_dark_B_fall_MC , total_contribution_pH_B_fall_MC , total_contribution_sun_B_fall_MC];

contribution_A_MC = [mean([contribution_all_A_winter_MC(:,1) , contribution_all_A_spring_MC(:,1) , contribution_all_A_summer_MC(:,1) , contribution_all_A_fall_MC(:,1)],2),...
    mean([contribution_all_A_winter_MC(:,2) , contribution_all_A_spring_MC(:,2) , contribution_all_A_summer_MC(:,2) , contribution_all_A_fall_MC(:,2)],2),...
    mean([contribution_all_A_winter_MC(:,3) , contribution_all_A_spring_MC(:,3) , contribution_all_A_summer_MC(:,3) , contribution_all_A_fall_MC(:,3)],2)];
contribution_B_MC = [mean([contribution_all_B_winter_MC(:,1) , contribution_all_B_spring_MC(:,1) , contribution_all_B_summer_MC(:,1) , contribution_all_B_fall_MC(:,1)],2),...
    mean([contribution_all_B_winter_MC(:,2) , contribution_all_B_spring_MC(:,2) , contribution_all_B_summer_MC(:,2) , contribution_all_B_fall_MC(:,2)],2),...
    mean([contribution_all_B_winter_MC(:,3) , contribution_all_B_spring_MC(:,3) , contribution_all_B_summer_MC(:,3) , contribution_all_B_fall_MC(:,3)],2)];

%% Figure 3

position_x_left = 0.10;
position_x_right = 0.55;
position_y_low = 0.05;
position_y_high = 0.53;
position_width = 0.40;
position_height = 0.38;

figure(101), clf, hold on
X = [0.75  2.5  4.25];


subplot(2,2,1), hold on
h1 = boxplot(contribution_all_A_winter_MC*100,'positions',X-0.25,...
    'colors','k','width',0.45,'PlotStyle','traditional',...
    'MedianStyle','line','OutlierSize',3,'whisker',0.7193);
h2 = boxplot(contribution_all_B_winter_MC*100,'positions',X+0.25,...
    'colors',[0.5 0.5 0.5],'width',0.45,'PlotStyle','traditional',...
    'MedianStyle','line','OutlierSize',3,'whisker',0.7193);
set(gca,'xticklabel',{'Dark decay' 'pH toxicity' 'Sunlight damage'},'FontSize',fs-2,...
    'Position',[position_x_left position_y_high position_width position_height])
xlim([0 5])
ylim([0 100])
title('Winter','FontSize',fs)
set(gca,'LineWidth',3)
set(h1,'LineWidth',2)
set(h2,'LineWidth',2)
yticks(0:20:100)
h = findobj(gca,'Tag','Box');
hLegend = legend(h([4,1]), {'HRAP A','HRAP B'},'FontSize',fs,'Location','North');

% scatter([1 2 3],[total_contribution_dark_A_winter,total_contribution_pH_A_winter,total_contribution_sun_A_winter]*100,100,'+r','LineWidth',1)


subplot(2,2,2), hold on
h1 = boxplot(contribution_all_A_spring_MC*100,'positions',X-0.25,...
    'colors','k','width',0.45,'PlotStyle','traditional',...
    'MedianStyle','line','OutlierSize',3,'whisker',0.7193);
h2 = boxplot(contribution_all_B_spring_MC*100,'positions',X+0.25,...
    'colors',[0.5 0.5 0.5],'width',0.45,'PlotStyle','traditional',...
    'MedianStyle','line','OutlierSize',3,'whisker',0.7193);
set(gca,'xticklabel',{'Dark decay' 'pH toxicity' 'Sunlight damage'},'FontSize',fs-2,...
    'Position',[position_x_right position_y_high position_width position_height])
xlim([0 5])
ylim([0 100])
title('Spring')
set(gca,'LineWidth',3)
set(h1,'LineWidth',2)
set(h2,'LineWidth',2)
yticks(0:20:100)
% scatter([1 2 3],[total_contribution_dark_A_spring,total_contribution_pH_A_spring,total_contribution_sun_A_spring]*100,100,'+r','LineWidth',1)


subplot(2,2,4), hold on
h1 = boxplot(contribution_all_A_summer_MC*100,'positions',X-0.25,...
    'colors','k','width',0.45,'PlotStyle','traditional',...
    'MedianStyle','line','OutlierSize',3,'whisker',0.7193);
h2 = boxplot(contribution_all_B_summer_MC*100,'positions',X+0.25,...
    'colors',[0.5 0.5 0.5],'width',0.45,'PlotStyle','traditional',...
    'MedianStyle','line','OutlierSize',3,'whisker',0.7193);
set(gca,'xticklabel',{'Dark decay' 'pH toxicity' 'Sunlight damage'},'FontSize',fs-2,...
    'Position',[position_x_right position_y_low position_width position_height])
xlim([0 5])
ylim([0 100])
title('Summer')
set(gca,'LineWidth',3)
set(h1,'LineWidth',2)
set(h2,'LineWidth',2)
yticks(0:20:100)
% scatter([1 2 3],[total_contribution_dark_A_summer,total_contribution_pH_A_summer,total_contribution_sun_A_summer]*100,100,'+r','LineWidth',1)


subplot(2,2,3), hold on
h1 = boxplot(contribution_all_A_fall_MC*100,'positions',X-0.25,...
    'colors','k','width',0.45,'PlotStyle','traditional',...
    'MedianStyle','line','OutlierSize',3,'whisker',0.7193);
h2 = boxplot(contribution_all_B_fall_MC*100,'positions',X+0.25,...
    'colors',[0.5 0.5 0.5],'width',0.45,'PlotStyle','traditional',...
    'MedianStyle','line','OutlierSize',3,'whisker',0.7193);
set(gca,'xticklabel',{'Dark decay' 'pH toxicity' 'Sunlight damage'},'FontSize',fs-2,...
    'Position',[position_x_left position_y_low position_width position_height])
xlim([0 5])
ylim([0 100])
title('Fall')
set(gca,'LineWidth',3)
set(h1,'LineWidth',2)
set(h2,'LineWidth',2)
yticks(0:20:100)
% scatter([1 2 3],[total_contribution_dark_A_fall,total_contribution_pH_A_fall,total_contribution_sun_A_fall]*100,100,'+r','LineWidth',1)


figure(101)
fig = gcf;
fig.Position = [50 50 900 700];
yll = ylabel('Mechanism relative contribution to total {\itE. coli} decay (%)','FontSize',fs + 1,'FontWeight','bold')
yll.Position(1) = yll.Position(1) - abs(yll.Position(1)*0.5)
yll.Position(2) = yll.Position(2) + abs(yll.Position(2)*1.25)

