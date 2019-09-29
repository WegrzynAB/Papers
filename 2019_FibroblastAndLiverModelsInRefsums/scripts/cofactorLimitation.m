%Programmer:   Agnieszka Wegrzyn
%Adapted from: Maria Suarez Diez, Rienk Rienksma Ironlimitation.m
%Script to calculate flux distributions on standard Recon when cofactor is limited  (FAD)

% 20190322: updated to replace artificial reaction with coupling constraint
% Agnieszka Wegrzyn
%% Initialize model
clear all;
initCobraToolbox(false);
modelOrig = 'Recon3'; %'Recon3' for Recon3D 'Recon2' for Recon2.2
File = 'E:/Dropbox/Cofactor Manuscript/scripts/models/modelR3D_X.mat'; 
load(File);
%% Read Genes encoding enzymes that use  cofactor (provide correct pathway to your local S1_Table)
if strcmp('Recon3', modelOrig)
        [~, ~, raw] = xlsread('E:/Dropbox/Cofactor Manuscript/Additional_files/S1_Table.xlsx','Flavoproteins','G2:G112');
        model = modelR3D_X;
    else
        [~, ~, raw] = xlsread('E:/Dropbox/Cofactor Manuscript/Additional_files/S1_Table.xlsx','Flavoproteins','J2:J112');
        model = modelR22_X;
end

raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,1);
flavoprot = cellVectors(:,1);
clearvars raw cellVectors;

%% Find the genes in the model associated with the cofactor genes
% clean up exchange, demand, or sink reaction for FAD and add a pseudo reaction converting FAD to cofactor
model = removeRxns(model, 'DM_fad[c]');
%model = removeRxns(model, 'sink_fad[c]');
model = removeRxns(model, 'EX_fad[c]');
model = removeRxns(model, 'EX_fad[e]');

%model= addReaction(model,'FADisCofactor', {'fad[c]', 'cofactor'}, [-1,1] ,0, 0, 1000);
%model.subSystems(end) = {''};
flavoprotMapped = flavoprot(ismember(flavoprot, model.genes));

% find the reactions in the model affected by these genes (consider boolean rules in the GPRs).
[~,~,constrRxnNames,~] = deleteModelGenes(model,flavoprotMapped);
% make the model irreversible regarding these reactions. 
% [model,flavoprotRxns] = partIrrev(model,constrRxnNames);
%% add coupling constraint
rxnList = model.rxns(findRxnIDs(model,constrRxnNames));
sol = optimizeCbModel(changeObjective(model, rxnList), 'max');
rxnC = 'FMNAT';
u = 0.000002;
c = (u-sol.f)/-5;

[modelCoupled] = coupleRxnList2Rxn(model, rxnList, rxnC, c, u);
[modelCoupled, removedRxnInd, keptRxnInd] = checkDuplicateRxn(modelCoupled, 'S', 1, 0);
%% make all flavoprotein-related reactions consume 'cofactor'
% Sendrow = model.S(end,:);
% for i = 1:1:length(flavoprotRxns)
%     if flavoprotRxns(i) == 1 
%         Sendrow(1,i) = -0.000002*0.1; %-0.000002;
%     end
% end
% 
% model.S(end,:) = Sendrow;

%% save results
if strcmp('Recon3', modelOrig)
        modelR3D_X_c = modelCoupled;
        save  'modelR3D_X_c.mat' modelR3D_X_c
        clearvars -except modelR3D modelR3D_X modelR3D_X_c flavoprotMapped;
    else
        modelR22_X_c = modelCoupled; 
        save  'modelR22_flavo_c.mat' modelR22_flavo_c 
        clearvars -except modelR22 modelR22_X modelR22_X_c flavoprotMapped;
end
%% test results
 CofSensTest_flavo;