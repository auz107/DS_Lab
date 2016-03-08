clear all

% Enable libsbml
addpath('/usr/local;/usr/local/lib;/usr/local/include');
addpath('/usr/local/builds/libSBML-5.10.2-Source/src/bindings/matlab/');

% Enable RAVEN
addpath('~/RAVEN')

%load('Pflur_iSB1139.mat')
model = importModel('iSB1139_orig.xml')
%model = importExcelModel('Pflur.xls',true)

% Find reactions with no associated reactants and products
counter = 0;
for i = 1:size(model.S,2)
    if length(find(model.S(:,i) ~= 0)) == 0
        fprintf('%s\n',model.rxns{i});
        counter = counter + 1;
    end
end

fprintf('\nThe total # of empty rxns = %d\n\n',counter)

% MIRXN_39: H+_input => H+[extracellular]
Hp_ext = find(strcmp(model.mets,'r718'))
