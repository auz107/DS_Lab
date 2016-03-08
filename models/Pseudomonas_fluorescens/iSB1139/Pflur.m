%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test code Partho
% 10/12/2014

clear all; close all; clc;

%model = importModel('Pflu.xml');
model=importExcelModel('Pflur_Partho.xls',true);
save('Pflur','model');

load('Pflur');

[exchangeRxns, exchangeRxnsIndexes] = getExchangeRxns(model,'UP');


 model.lb(exchangeRxnsIndexes) = -1000 %% change -> here the biomass changes
 model.ub(exchangeRxnsIndexes) = 1000   %% change -> here the biomass changes

solveLP(model,'max')
