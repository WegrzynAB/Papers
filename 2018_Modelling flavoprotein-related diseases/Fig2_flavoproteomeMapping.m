%% Initialize model
initCobraToolbox;
load('modelR22.mat')
load('modelR22_FAD.mat')
%% Read Genes encoding enzymes that use  cofactor (provide correct pathway to your local S1_Table)
[~, ~, raw] = xlsread('S1_Table.xlsx','Flavoproteins','J2:J112');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,1);
flavoprot = cellVectors(:,1);
clearvars raw cellVectors;
[~, ~, raw] = xlsread('S1_Table.xlsx','Flavoproteins','D2:D112');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,1);
flavoprot_Names = cellVectors(:,1);
clearvars raw cellVectors;
%% Find rxns associated with flavoproteins
genesR22 = flavoprot(ismember(flavoprot, modelR22.genes));
genesR22_F = flavoprot(ismember(flavoprot, modelR22_FAD.genes));
unique_genes = unique([genesR22; genesR22_F]);
unique_names = cell(length(unique_genes),1);
for i=1:length(unique_genes)
    unique_names(i) = flavoprot_Names(ismember(flavoprot, unique_genes(i)));
end
countsR22 = zeros(length(unique_genes),1);
countsR22_F = zeros(length(unique_genes),1);

for i=1:length(unique_genes)
    if ismember(unique_genes(i), modelR22.genes)
        [~,solR22] = findRxnsFromGenes(modelR22, unique_genes(i),0,1);
        countsR22(i) = length(solR22(:,1));
    end
    if ismember(unique_genes(i), modelR22_FAD.genes)
        [~,solR22_F] = findRxnsFromGenes(modelR22_FAD, unique_genes(i),0,1);
        countsR22_F(i) = length(solR22_F(:,1));
    end
end
%% results table for figure 2A
results_A = table(unique_genes, unique_names, countsR22, countsR22_F);

%% find in which subsystems flavoporoteins play an important role

unique_subSystems = unique([modelR22.subSystems;modelR22_FAD.subSystems]);
totalR22 = zeros(length(unique_subSystems),1);
totalR22_F = zeros(length(unique_subSystems),1);
flavoR22 = zeros(length(unique_subSystems),1);
flavoR22_F = zeros(length(unique_subSystems),1);
ratioR22 = zeros(length(unique_subSystems),1);
ratioR22_F = zeros(length(unique_subSystems),1);
[~,listR22] = findRxnsFromGenes(modelR22, flavoprot,0,1);
flavoRxnsR22 = unique(listR22);
[~,listR22_F] = findRxnsFromGenes(modelR22_FAD, flavoprot,0,1);
flavoRxnsR22_F = unique(listR22_F);

for i=1:length(unique_subSystems)
    solR22 = findRxnsFromSubSystem(modelR22, unique_subSystems(i));
    totalR22(i,1) = length(solR22(:,1));
    subSolR22 = solR22(ismember(solR22(:,1),flavoRxnsR22(:,1)));
    flavoR22(i,1) = length(subSolR22);
    ratioR22(i,1) = flavoR22(i,1)*100/totalR22(i,1);
    solR22_F = findRxnsFromSubSystem(modelR22_FAD, unique_subSystems(i));
    totalR22_F(i,1) = length(solR22_F(:,1));
    subSolR22_F = solR22_F(ismember(solR22_F(:,1),flavoRxnsR22_F(:,1)));
    flavoR22_F(i,1) = length(subSolR22_F);
    ratioR22_F(i,1) = flavoR22_F(i,1)*100/totalR22_F(i,1);
end

%% results table for figure 2B
results_B = table(unique_subSystems, flavoR22, totalR22, flavoR22_F, totalR22_F, ratioR22, ratioR22_F);

