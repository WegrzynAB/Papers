%Programmer:   Agnieszka Wegrzyn
%Adapted from: Maria Suarez Diez, Rienk Rienksma Ironlimitation.m
%Last updated: January 2017
%Script to calculate flux distributions on standard Recon when cofactor 
%is limited  (FAD)
%% Initialize model
%clear all;
%initCobraToolbox;
%load 'Recon22model_FAD.mat'
model = modelR22_FAD;

%% Read Genes encoding enzymes that use  cofactor (provide pathway to Supplementary Table 1 if needed)
[~, ~, raw] = xlsread('/Users/Dave/Dropbox/Cofactor Manuscript/Additional_files/S1_Table.xlsx','Flavoproteins','J2:J112');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,1);
flavoprot = cellVectors(:,1);
clearvars raw cellVectors;

%% Find the genes in the model associated with the cofactor genes
% add a pseudo reaction that measures the flux through the network
model = removeRxns(model, 'DM_fad[c]');
model = removeRxns(model, 'EX_fad[c]');
model= addReaction(model,'FADisCofactor', {'fad[c]', 'cofactor'}, [-1,1] ,[0], 0, 1000);

flavoprotMapped = flavoprot(ismember(flavoprot, model.genes));

% find the reactions in the model affected by these genes (consider
% boolean rules in the GPRs. 
[~,~,constrRxnNames,~]=deleteModelGenes(model,flavoprotMapped);
% make the model irreversible regarding these reactions. 
[model,flavoprotRxns] = partIrrev(model,constrRxnNames);
%% make all flavoprotein-related reactions consume 'cofactor'
Sendrow = model.S(end,:);
for i = 1:1:length(flavoprotRxns)
    if flavoprotRxns(i) == 1 
        Sendrow(1,i) = -0.000002;
    end
end

model.S(end,:) = Sendrow;
modelR22_flavo = model;
save  'modelR22_flavo.mat' modelR22_flavo
save  'modelR22_FAD.mat' modelR22_FAD