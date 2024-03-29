%Programmer:   Agnieszka Wegrzyn
%Adapted from: Maria Suarez Diez, Rienk Rienksma Ironlimitation.m
%Script to calculate flux distributions on standard Recon when cofactor is limited  (FAD)

%% Initialize model
initCobraToolbox;
File = 'modelR22_FAD.mat'; %'modelR3D_FAD.mat' for Recon3D 'modelR22_FAD.mat' for Recon2.2
load(File);
%% Read Genes encoding enzymes that use  cofactor (provide correct pathway to your local S1_Table)
if strcmp('modelR3D_FAD.mat',File)
        [~, ~, raw] = xlsread('S1_Table.xlsx','Flavoproteins','G2:G112');
        model = modelR3D_FAD;
    else
        [~, ~, raw] = xlsread('S1_Table.xlsx','Flavoproteins','J2:J112');
        model = modelR22_FAD;
end

raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,1);
flavoprot = cellVectors(:,1);
clearvars raw cellVectors;

%% Find the genes in the model associated with the cofactor genes
% clean up exchange, demand, or sink reaction for FAD and add a pseudo reaction converting FAD to cofactor
model = removeRxns(model, 'DM_fad[c]');
model = removeRxns(model, 'sink_fad[c]');
model = removeRxns(model, 'EX_fad[c]');
model = removeRxns(model, 'EX_fad[e]');

model= addReaction(model,'FADisCofactor', {'fad[c]', 'cofactor'}, [-1,1] ,0, 0, 1000);
%model.subSystems(end) = {''};
flavoprotMapped = flavoprot(ismember(flavoprot, model.genes));

% find the reactions in the model affected by these genes (consider boolean rules in the GPRs).
[~,~,constrRxnNames,~] = deleteModelGenes(model,flavoprotMapped);
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

%% save results
if strcmp('modelR3D_FAD.mat',File)
        modelR3D_flavo = model;
        save  'modelR3D_flavo.mat' modelR3D_flavo
        clearvars -except modelR3D modelR3D_FAD modelR3D_flavo flavoprotMapped;
    else
        modelR22_flavo = model; 
        save  'modelR22_flavo.mat' modelR22_flavo 
        clearvars -except modelR22 modelR22_FAD modelR22_flavo flavoprotMapped;
end
