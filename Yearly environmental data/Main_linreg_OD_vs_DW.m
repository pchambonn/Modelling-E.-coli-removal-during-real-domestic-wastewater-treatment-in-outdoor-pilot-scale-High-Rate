%% linreg_OD_vs_DW
% This script aims at performing the linear regression between optical 
% density and dryweight using data from pilot scale HRAP in order to
% determine the regression coefficient and accuracy on their determination
% for the publication 'Modelling E. coli removal during real domestic
% wastewater treatment in outdoor pilot scale High Rate Algal Ponds'.

% This script contributes to the results discussed in S3

%%
% Data import
A = xlsread('./Yearly environmental data/Optical_density_vs_TSS.xlsx')

TSS = A(:,1);
OD = A(:,4);

% Linear Regression
linearModel = fitlm(TSS,OD)
linearModel.Coefficients.SE*tinv(0.025,49)