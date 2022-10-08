%% Main_yearly_environmental_data_figure

% This script aims at loading the pre-existing data variables constructed
% with Main_yearly_environmental_data_format.m and at drawing the figures to
% be included in the publication PREDICTION OF E. COLI DECAY MECHANISMS
% DURING WASTEWATER TREATMENT IN HIGH RATE ALGAL PONDS relative to this
% data set.

% This script is hence used to draw the figure Figure S9.1 and carry out
% the statistical analysis used in the discussion of S9.


%% Data import
clear all
load('./Yearly environmental data/yearly_environmental_data_formatted.mat')
fs = 12;

%% Environmental data boxplots (Figure S9.1)

figure(1), clf, hold on

fig = gcf;
fig.Position = [30 40 900 700];

X = [1 2.5 4 5.5 7];

subplot(2,2,1) % sunlight
h1 = boxplot(Hs_plot,'positions',X,'colors','k','width',0.45,'PlotStyle','traditional','MedianStyle','line','OutlierSize',3,'whisker',0.7193);
ylabel('Sunlight intensity (W.m^{-2})')
xticklabels([])
ylim([0 1200])
xlim([min(X) - 1 max(X) + 1]);
ax = gca;
ax.Position = [0.09 0.55 0.39 0.39];
xticks(ax,X)
ax.FontSize = fs; ax.FontWeight = 'bold';
ax2 = axes('Position',get(ax,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','Box','off');
xlim(ax.XLim); xticks(ax.XTick); set(ax2,'YTick',[],'YTickLabel',[])  
xticklabels({num2str(sum(abs(1 - isnan(Hs_plot(:,1))))),...
             num2str(sum(abs(1 - isnan(Hs_plot(:,2))))),...
             num2str(sum(abs(1 - isnan(Hs_plot(:,3))))),...
    	     num2str(sum(abs(1 - isnan(Hs_plot(:,4))))),...
             num2str(sum(abs(1 - isnan(Hs_plot(:,5)))))})
ax2.FontSize = fs - 4;
set(ax,'LineWidth',3)
set(h1,'LineWidth',2)


subplot(2,2,2), hold on % sunlight
h1 = boxplot(T_A_plot,'positions',X - 0.3,'colors','k','width',0.45,'PlotStyle','traditional','MedianStyle','line','OutlierSize',3,'whisker',0.7193);
h2 = boxplot(T_B_plot,'positions',X + 0.3,'colors',[0.5 0.5 0.5],'width',0.45,'PlotStyle','traditional','MedianStyle','line','OutlierSize',3,'whisker',0.7193);
ylabel('Broth temperature (°C)')
xticklabels([])
xlim([min(X) - 1 max(X) + 1])
ylim([0 30])
ax = gca;
ax.Position = [0.59 0.55 0.39 0.39];
xticks(ax,X)
ax.FontSize = fs; ax.FontWeight = 'bold';
ax2 = axes('Position',get(ax,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','Box','off');
xlim(ax.XLim); set(ax2,'YTick',[],'YTickLabel',[])
xticks(sort(cat(2,X - 0.4,X,X + 0.4)))
xticklabels({num2str(sum(abs(1 - isnan(T_A_plot(:,1))))),'/',num2str(sum(abs(1 - isnan(T_B_plot(:,1))))),...
             num2str(sum(abs(1 - isnan(T_A_plot(:,2))))),'/',num2str(sum(abs(1 - isnan(T_B_plot(:,2))))),...
             num2str(sum(abs(1 - isnan(T_A_plot(:,3))))),'/',num2str(sum(abs(1 - isnan(T_B_plot(:,3))))),...
             num2str(sum(abs(1 - isnan(T_A_plot(:,4))))),'/',num2str(sum(abs(1 - isnan(T_B_plot(:,4))))),...
             num2str(sum(abs(1 - isnan(T_A_plot(:,5))))),'/',num2str(sum(abs(1 - isnan(T_B_plot(:,5)))))})
ax2.FontSize = fs - 4;
set(ax,'LineWidth',3)
set(h1,'LineWidth',2)
set(h2,'LineWidth',2)


subplot(2,2,3), hold on % sunlight
h1 = boxplot(pH_A_plot,'positions',X - 0.3,'colors','k','width',0.45,'PlotStyle','traditional','MedianStyle','line','OutlierSize',3,'whisker',0.7193);
h2 = boxplot(pH_B_plot,'positions',X + 0.3,'colors',[0.5 0.5 0.5],'width',0.45,'PlotStyle','traditional','MedianStyle','line','OutlierSize',3,'whisker',0.7193);
ylabel('Broth pH')
ax = gca;
ax.Position = [0.09 0.05 0.39 0.39];
xticks(ax,X)
xticklabels({'Overall' , 'Winter', 'Spring' , 'Summer' , 'Fall'})
xlim([min(X) - 1 max(X) + 1])
ylim([6.5 11.5])
ax.FontSize = fs; ax.FontWeight = 'bold';
ax2 = axes('Position',get(ax,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','Box','off');
xlim(ax.XLim);set(ax2,'YTick',[],'YTickLabel',[])  
xticks(sort(cat(2,X - 0.4,X,X + 0.4)))
xticklabels({num2str(sum(abs(1 - isnan(pH_A_plot(:,1))))),'/',num2str(sum(abs(1 - isnan(pH_B_plot(:,1))))),...
             num2str(sum(abs(1 - isnan(pH_A_plot(:,2))))),'/',num2str(sum(abs(1 - isnan(pH_B_plot(:,2))))),...
             num2str(sum(abs(1 - isnan(pH_A_plot(:,3))))),'/',num2str(sum(abs(1 - isnan(pH_B_plot(:,3))))),...
             num2str(sum(abs(1 - isnan(pH_A_plot(:,4))))),'/',num2str(sum(abs(1 - isnan(pH_B_plot(:,4))))),...
             num2str(sum(abs(1 - isnan(pH_A_plot(:,5))))),'/',num2str(sum(abs(1 - isnan(pH_B_plot(:,5)))))})
ax2.FontSize = fs - 4;
set(ax,'LineWidth',3)
set(h1,'LineWidth',2)
set(h2,'LineWidth',2)


subplot(2,2,4), hold on % sunlight
h1 = boxplot(sigma_A_plot,'positions',X - 0.3,'colors','k','width',0.45,'PlotStyle','traditional','MedianStyle','line','OutlierSize',3,'whisker',0.7193);
h2 = boxplot(sigma_B_plot,'positions',X + 0.3,'colors',[0.5 0.5 0.5],'width',0.45,'PlotStyle','traditional','MedianStyle','line','OutlierSize',3,'whisker',0.7193);
ylabel({'Broth light attenuation' ;  'coefficient (m^{-1})'})
xticklabels({'Overall' , 'Winter', 'Spring' , 'Summer' , 'Fall'})
ax = gca;
ax.Position = [0.59 0.05 0.39 0.39];
xticks(ax,X)
xlim([min(X) - 1 max(X) + 1])
ylim([0 200])
ax.FontSize = fs; ax.FontWeight = 'bold';
ax2 = axes('Position',get(ax,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','Box','off');
xlim(ax.XLim); set(ax2,'YTick',[],'YTickLabel',[])  
xticks(sort(cat(2,X - 0.4,X,X + 0.4)))
xticklabels({num2str(sum(abs(1 - isnan(sigma_A_plot(:,1))))),'/',num2str(sum(abs(1 - isnan(sigma_B_plot(:,1))))),...
             num2str(sum(abs(1 - isnan(sigma_A_plot(:,2))))),'/',num2str(sum(abs(1 - isnan(sigma_B_plot(:,2))))),...
             num2str(sum(abs(1 - isnan(sigma_A_plot(:,3))))),'/',num2str(sum(abs(1 - isnan(sigma_B_plot(:,3))))),...
             num2str(sum(abs(1 - isnan(sigma_A_plot(:,4))))),'/',num2str(sum(abs(1 - isnan(sigma_B_plot(:,4))))),...
             num2str(sum(abs(1 - isnan(sigma_A_plot(:,5))))),'/',num2str(sum(abs(1 - isnan(sigma_B_plot(:,5)))))})
ax2.FontSize = fs - 4;
set(ax,'LineWidth',3)
set(h1,'LineWidth',2)
set(h2,'LineWidth',2)


h = findobj(ax,'Tag','Box');
hLegend = legend(h([6,1]), {'HRAP A','HRAP B'},'FontSize',fs);

%% Statistical tests for similar conditions
% The statistical analysis detailed in Supplementary Information (S9) is
% performed here. The objective was to evaluate if the environmental
% conditions experienced in both ponds were comparable to help explain
% potential discrepancies in the decay right predicted. Environmental
% variables (e.g. broth temperature) were comapred in terms of means,
% and standard deviations.
% As explained in the manuscript, equality of means was tested using
% two-samples t-tests, and equality of variance was tested using Bartlett
% test. Additional one tailed two-samples Kolmogorov Smirnov tests were
% carried out to evidence if values experienced in a pilot HRAP was higher
% or lower than in the other.
        

% Temperature
T_A = data_A(:,3);
T_B = data_B(:,3);

[h_ttest_T,p_ttest_T,ci_ttest_T,stats_ttest_T] = ttest2(T_A,T_B,'tail','both','Vartype','unequal')

T_A_bartlett = [T_A;NaN(max(length(T_B) - length(T_A),0),1)];
T_B_bartlett = [T_B;NaN(max(length(T_A) - length(T_B),0),1)];
T_bartlett = [T_A_bartlett,T_B_bartlett];
[p_bartlett_T,stats_bartlett_T] = vartestn(T_bartlett)

[h_ks_T,p_ks_T,stats_ks_T] = kstest2(T_A,T_B,'tail','larger')


% pH
pH_A = data_A(:,2);
pH_B = data_B(:,2);

[h_ttest_pH,p_ttest_pH,ci_ttest_pH,stats_ttest_pH] = ttest2(pH_A,pH_B,'tail','both','Vartype','unequal')

pH_A_bartlett = [pH_A;NaN(max(length(pH_B) - length(pH_A),0),1)];
pH_B_bartlett = [pH_B;NaN(max(length(pH_A) - length(pH_B),0),1)];
pH_bartlett = [pH_A_bartlett,pH_B_bartlett];
[p_bartlett_pH,stats_bartlett_pH] = vartestn(pH_bartlett)

[h_ks_pH,p_ks_pH,stats_ks_pH] = kstest2(pH_A,pH_B,'tail','larger')


% sigma
data_sigma = xlsread('./Yearly environmental data/TSS_light attenuation.xlsx');
sigma_A_meas = data_sigma(:,3);
sigma_A_meas = sigma_A_meas(isnan(sigma_A_meas) == 0);
sigma_B_meas = data_sigma(:,6);
sigma_B_meas = sigma_B_meas(isnan(sigma_B_meas) == 0);

[h_ttest_sigma,p_ttest_sigma,ci_ttest_sigma,stats_ttest_sigma] = ttest2(sigma_A_meas,sigma_B_meas,'tail','both','Vartype','unequal')

sigma_A_bartlett = [sigma_A_meas(:,1);NaN(max(length(sigma_B_meas) - length(sigma_A_meas),0),1)];
sigma_B_bartlett = [sigma_B_meas(:,1);NaN(max(length(sigma_A_meas) - length(sigma_B_meas),0),1)];
sigma_bartlett = [sigma_A_bartlett,sigma_B_bartlett];
[p_bartlett_sigma,stats_bartlett_sigma] = vartestn(sigma_bartlett)

[h_ks_sigma,p_ks_sigma,stats_ks_sigma] = kstest2(sigma_A_meas,sigma_B_meas,'tail','larger')




