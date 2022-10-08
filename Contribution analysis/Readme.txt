This folder contains the formatted data files and scripts used to analyze mechanism contribution to E. coli decay in the publication 'Modelling E. coli removal during real domestic wastewater treatment in outdoor pilot scale High Rate Algal Ponds'.

Matlab data sets:
* yearlyDecayContribution.mat:
	Data set of values for mechanism contribution to overall E. coli decay as predicted from the yearly environmental data.
* yearlyDecayContribution_MonteCarlo.mat:
	Data set of values for mechanism contribution to overall E. coli decay as predicted from the yearly environmental data. This data set includes 2000 data sets as generated from the 2000 sets of fitting parameters obtained from Monte Carlo analysis to evaluate uncertainties on predictions.
* yearlyDecayValue.mat:
	Data set of values for mechanism contribution to overall E. coli decay as predicted from the yearly environmental data. 
/!\  the data sets were not included to the GitHub repository. These data sets should therefore be generated first from the available scripts before running the scripts needing this data


Matla scripts:
* Main_plot_seasonalZoomedDecay.m:
	This script enables to plot Figure 4 and Figure S10.1 of the publication 'Modelling E. coli removal during real domestic wastewater treatment in outdoor pilot scale High Rate Algal Ponds'.
* Main_plot_yearlySeasonalContribution.m: 
	This script enables to plot Figure 3 and Figure S9.2 of the publication 'Modelling E. coli removal during real domestic wastewater treatment in outdoor pilot scale High Rate Algal Ponds'.
* Main_yearlyMechanismContribution_bestFit.m: 
	This script aims at computing the decay CONTRIBUTION of each mechanism to overall E. coli decay in a theoretical HRAP which environmental conditions were as measured i HRAP A and B. These vectors serve as the basis for plotting Figure S9.2 of the publication 'Modelling E. coli removal during real domestic wastewater treatment in outdoor pilot scale High Rate Algal Ponds' (as performed in the script Main_plot_yearlySeasonalContribution). 
* Main_yearlyMechanismContribution_monteCarloUncertainty.m:
	THis script aims at computing the contribution of each mechanism to overall E. coli decay in a theoretical HRAP which environmental conditions were as measured i HRAP A and B including variations for every Monte Carlo set of value for fitting parameters. These vectors serve as the basis for plotting Figure S9.2 of the publication 'Modelling E. coli removal during real domestic wastewater treatment in outdoor pilot scale High Rate Algal Ponds' (as performed in the script Main_plot_yearlySeasonalContribution). 
* Main_yearlyMechanismValue_bestFit.m:
	THis script aims at computing the decay VALUE of each mechanism to overall E. coli decay in a theoretical HRAP which environmental conditions were as measured i HRAP A and B. These vectors serve as the basis for plotting Figure 4 and Figure S10.1 of the publication 'Modelling E. coli removal during real domestic wastewater treatment in outdoor pilot scale High Rate Algal Ponds' (as performed in the script Main_plot_seasonalZoomedDecay). 
