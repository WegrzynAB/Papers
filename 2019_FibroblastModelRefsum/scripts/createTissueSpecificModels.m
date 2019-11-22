%% create tissue specific models (fibroblast and liver)
%initCobraToolbox
load('E:\Dropbox\Refsum\manuscript\models\modelR3D_FAD_X.mat')
importMediaConstraints
importFibroblastConstraints
model = changeObjective(modelR3D_FAD_X, 'EX_phyt[e]');
model = addExchangeRxn(model, 'zn2[e]');
model = addExchangeRxn(model, 'ca2[e]');
model = addExchangeRxn(model, 'cl[e]');
modelExRxns = model.rxns(findExcRxns(model,1));
activeImport = modelExRxns(ismember(modelExRxns, findRxnsFromMets(model, ModelMetIDs)));
%inactiveImport = modelExRxns(~ismember(modelExRxns, findRxnsFromMets(model, ModelMetIDs)));
model = blockAllImports(model);
simple_media = {'EX_ca2[e]'; 'EX_cl[e]'; 'EX_fe2[e]'; 'EX_fe3[e]';...
    'EX_h[e]'; 'EX_h2o[e]'; 'EX_k[e]'; 'EX_na1[e]'; 'EX_nh4[e]';...
    'EX_so4[e]'; 'EX_pi[e]'; 'EX_o2[e]'};
model = changeRxnBounds(model, simple_media, -1000, 'l');
model = changeRxnBounds(model, activeImport, lb, 'l');
decision = logical(DECISION);

for i = 1:length(EntrezID)
    gene = num2str(EntrezID(i));
    genes{i,1} = {gene};
end
activeGenes = genes(decision);
modelGenes = model.genes;
inactiveGenes = genes(~decision);
[~, ~, inIDs] = findRxnsFromEntrezIDs(model, inactiveGenes);
[activeRxnsList, geneIdRxns, IDs] = findRxnsFromEntrezIDs(model, activeGenes);
activeRxnsList = [activeRxnsList; ActiveRxns];
activeRxnsList = [activeRxnsList; activeImport];
activeRxnsList = unique(activeRxnsList);

model = changeRxnBounds(model, 'DM_atp_c_', 0.1, 'l');
sol = optimizeCbModel(changeObjective(model, 'biomass_maintenance'));
model = changeRxnBounds(model, 'biomass_maintenance', 0.01*sol.f, 'l'); %fix lb as maintenance
model = changeRxnBounds(model, 'PHYHx', 0.1, 'l');
model = changeRxnBounds(model, 'EX_3MAA[e]', 0.1, 'l');
model = changeRxnBounds(model, 'EX_but[e]', -0.1, 'u');
model = changeRxnBounds(model, 'EX_caproic[e]', -10, 'l');
model = changeRxnBounds(model, 'EX_caproic[e]', -0.1, 'u');
model = changeRxnBounds(model, 'EX_dca[e]', -0.1, 'b');
model = changeRxnBounds(model, 'EX_hexc[e]', -0.1, 'b');
model = changeRxnBounds(model, 'EX_lgnc[e]', -0.1, 'b');
model = changeRxnBounds(model, 'EX_docosac[e]', -0.1, 'b');
model = changeRxnBounds(model, 'EX_arach[e]', -0.1, 'b');
model = changeRxnBounds(model, 'EX_ocdca[e]', -0.1, 'b');
model = changeRxnBounds(model, 'EX_hdca[e]', -0.1, 'b');
model = changeRxnBounds(model, 'EX_ttdca[e]', -0.1, 'b');
model = changeRxnBounds(model, 'EX_ddca[e]', -0.1, 'b');
epsilon = 1e-4;
[A] = fastcc(model, epsilon, 2);
consistRxnsModel = model.rxns(A);
inconsistRxns = setdiff(model.rxns,consistRxnsModel);
modelConsistent = removeRxns(model,inconsistRxns);
[A] = fastcc(modelConsistent, epsilon, 2);
%modelConsistent = addExchangeRxn(modelConsistent, 'ps_hs[e]');
%define core reactions from model
coreConsistentRxns = find(ismember(modelConsistent.rxns,activeRxnsList));
coreConsistentRxns = sort(coreConsistentRxns);
% solve problem by finding the most compact subnetwork containing all core reactions
A2 = fastcore(modelConsistent, coreConsistentRxns, epsilon, 2);
%remove the reactions not in the subnetwork generated by fastcore
ind = findRxnIDs(modelConsistent, A2.rxns);
cellModelRxns = modelConsistent.rxns(ind);
notCellRxns = setdiff(modelConsistent.rxns,cellModelRxns);
modelSpecific = removeRxns(modelConsistent,notCellRxns);
metIDs = [];
for i=1:length(ModelMetIDs)
    metID = findMetIDs(model, ModelMetIDs(i));
    metIDs = [metIDs; metID];
end
metNames = model.mets(metIDs);
modelSpecific = addExchangeRxn(modelSpecific, metNames);
modelSpecific = addExchangeRxn(modelSpecific, 'h2o[e]');
modelSpecific = addExchangeRxn(modelSpecific, 'caproic[e]');
modelSpecific = addExchangeRxn(modelSpecific, 'ps_hs[e]');
[A] = fastcc(modelSpecific, epsilon, 2);
consistRxnsModel = modelSpecific.rxns(A);
inconsistRxns = setdiff(modelSpecific.rxns,consistRxnsModel);
modelSpecific = removeRxns(modelSpecific,inconsistRxns);

save('model_temp.mat', 'modelSpecific')

inactiveModelGenes = modelSpecific.genes(inIDs);
[~, ~, constrRxnNames, ~] = deleteModelGenes(modelSpecific, inactiveModelGenes);
inactiveRxnsList = constrRxnNames(~ismember(constrRxnNames, activeRxnsList));
modelDel = removeRxns(modelSpecific,inactiveRxnsList);

[A] = fastcc(modelDel, epsilon, 2);
consistRxnsModel = modelDel.rxns(A);
inconsistRxns = setdiff(modelDel.rxns,consistRxnsModel);
modelSpecific = removeRxns(modelSpecific,inconsistRxns);

save('model_temp2.mat', 'modelSpecific')

modelSpecific = changeRxnBounds(modelSpecific, 'DM_atp_c_', 0, 'l');
sol = optimizeCbModel(changeObjective(modelSpecific, 'biomass_maintenance'));
modelSpecific = changeRxnBounds(modelSpecific, 'biomass_maintenance', 0.01*sol.f, 'l'); %fix lb as maintenance
modelSpecific = changeRxnBounds(modelSpecific, 'PHYHx', 0, 'l');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_3MAA[e]', 0, 'l');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_but[e]', 1000, 'u');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_caproic[e]', 1000, 'u');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_dca[e]', 1000, 'u');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_dca[e]', 0, 'l');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_hexc[e]', 0, 'l');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_lgnc[e]', 0, 'l');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_docosac[e]', 0, 'l');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_arach[e]', 0, 'l');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_ocdca[e]',0, 'l');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_hdca[e]', 0, 'l');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_ttdca[e]', 0, 'l');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_ddca[e]', 0, 'l');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_hexc[e]', 1000, 'u');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_lgnc[e]', 1000, 'u');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_docosac[e]', 1000, 'u');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_arach[e]', 1000, 'u');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_ocdca[e]', 1000, 'u');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_hdca[e]', 1000, 'u');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_ttdca[e]', 1000, 'u');
modelSpecific = changeRxnBounds(modelSpecific, 'EX_ddca[e]', 1000, 'u');
modelSpecific = checkDuplicateRxn(modelSpecific, 'S', 1, 0);
modelSpecific = changeRxnBounds(modelSpecific, 'HMR_0180', -1000, 'l');

[A] = fastcc(modelSpecific, epsilon, 2);
consistRxnsModel = modelSpecific.rxns(A);
inconsistRxns = setdiff(modelSpecific.rxns,consistRxnsModel);
modelSpecific = removeRxns(modelSpecific,inconsistRxns);
modelSpecific = removeUnusedGenes(modelSpecific);

modelSpecific = addExchangeRxn(modelSpecific, {'co2[e]', 'nh4[e]'});
fibroblast = modelSpecific;

save('fibroblast.mat', 'fibroblast')
%% test model for leaks, metabolic functions, etc.
model= changeRxnBounds(fibroblast, 'biomass_maintenance', 0, 'l');
tutorial_modelSanityChecks_2
%% ATP yield
fibroblast = changeRxnBounds(fibroblast, 'EX_glc_D[e]', 0, 'u');
fibroblast = changeRxnBounds(fibroblast, 'biomass_maintenance', 0, 'l');
fibroblast = changeRxnBounds(fibroblast, 'DM_atp_c_', 0, 'l');
fibroblast_C = changeRxnBounds(fibroblast, 'CYP450phyt', 20.2176, 'u');
fibroblast_C = changeRxnBounds(fibroblast_C, 'PHYHx', 48.7656, 'u');
fibroblast_R = changeRxnBounds(fibroblast_C, 'PHYHx', 0, 'b');

[modelR3D_FAD_ATP] = maxFluxesB2(modelR3D_FAD,0);
[modelR3D_X_c_ATP] = maxFluxesB2(modelR3D_X_c,0);
[fibroblastC_ATP] = maxFluxesB2(fibroblast_C,0);
[fibroblastR_ATP] = maxFluxesB2(fibroblast_R,0);

summarised_table = table(modelR3D_FAD_ATP.carbon_source_t, modelR3D_FAD_ATP.flux_t, modelR3D_X_c_ATP.flux_t, fibroblastC_ATP.flux_t, fibroblastR_ATP.flux_t);
summarised_table.Properties.VariableNames = {'Carbon_source' 'modelR3D_FAD' 'modelR3D_X_c' 'fibroblastC' 'fibroblastR'};
writetable(summarised_table,'summarised_table_ATP_ext5.txt');

%% sampling of the fibroblast_C, fibroblast_R, fibroblast_C_phyt, fibroblast_R_phyt models
sol = optimizeCbModel(changeObjective(fibroblast, 'biomass_maintenance'));
fibroblast = changeRxnBounds(fibroblast, 'biomass_maintenance', 0.01*sol.f, 'l'); %fix lb as maintenance
fibroblast_C = changeRxnBounds(fibroblast, 'CYP450phyt', 20.2176, 'u');
fibroblast_C = changeRxnBounds(fibroblast_C, 'PHYHx', 48.7656, 'u');
options.nFiles = 10;
options.nPointsPerFile = 5000;
options.nPointsReturned = 1e4;
options.nStepsPerPoint = 500;
options.nFilesSkipped = 0;
options.nPointsPerFileLoaded = 1000;
options.maxTime = 16 * 3600;
tic
warmupPts= createHRWarmup(fibroblast_C);
ACHRSampler(fibroblast_C, warmupPts, 'control_fibro', options.nFiles, options.nPointsPerFile, options.nStepsPerPoint, [], [], options.maxTime, 0);
samples_C = loadSamples('control_fibro', 10, 1000, 0);
samples_C = samples_C(:, round(linspace(1, size(samples_C, 2), min([1e4, size(samples_C, 2)]))));
toc
save('fibroblast_C.mat', 'fibroblast_C', 'samples_C')
fibroblast_R = changeRxnBounds(fibroblast_C, 'PHYHx', 0, 'b');
tic
warmupPts = createHRWarmup(fibroblast_R);
ACHRSampler(fibroblast_R, warmupPts, 'refsum_fibro', options.nFiles, options.nPointsPerFile, options.nStepsPerPoint, [], [], options.maxTime, 0);
samples_R = loadSamples('refsum_fibro',  10, 1000, 0);
samples_R = samples_R(:, round(linspace(1, size(samples_R, 2), min([options.nPointsReturned, size(samples_R, 2)]))));
toc
save('fibroblast_R.mat', 'fibroblast_R', 'samples_R')
fibroblast_C_phyt = changeRxnBounds(fibroblast, 'EX_phyt[e]', -0.1, 'u');
fibroblast_R_phyt = changeRxnBounds(fibroblast_R, 'EX_phyt[e]', -0.1, 'u');
tic
warmupPts= createHRWarmup(fibroblast_C_phyt);
ACHRSampler(fibroblast_C_phyt, warmupPts, 'control_fibro_phyt', options.nFiles, options.nPointsPerFile, options.nStepsPerPoint, [], [], options.maxTime, 0);
samples_C_phyt = loadSamples('control_fibro_phyt', 10, 1000, 0);
samples_C_phyt = samples_C_phyt(:, round(linspace(1, size(samples_C_phyt, 2), min([options.nPointsReturned, size(samples_C_phyt, 2)]))));
toc
save('fibroblast_C_phyt.mat', 'fibroblast_C_phyt', 'samples_C_phyt')
fibroblast_R_phyt = changeRxnBounds(fibroblast_R_phyt, 'PHYHx', 0, 'b');
tic
warmupPts = createHRWarmup(fibroblast_R_phyt);
ACHRSampler(fibroblast_R_phyt, warmupPts, 'refsum_fibro_phyt', options.nFiles, options.nPointsPerFile, options.nStepsPerPoint, [], [], options.maxTime, 0);
samples_R_phyt = loadSamples('refsum_fibro_phyt', 10, 1000, 0);
samples_R_phyt = samples_R_phyt(:, round(linspace(1, size(samples_R_phyt, 2), min([options.nPointsReturned, size(samples_R_phyt, 2)]))));
toc
save('fibroblast_R_phyt.mat', 'fibroblast_R_phyt', 'samples_R_phyt')

csvwrite('sampledPoints_C.csv', samples_C)
csvwrite('sampledPoints_R.csv', samples_R)
csvwrite('sampledPoints_C_phyt.csv', samples_C_phyt)
csvwrite('sampledPoints_R_phyt.csv', samples_R_phyt)
