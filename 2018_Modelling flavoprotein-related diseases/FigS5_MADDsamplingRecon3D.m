load('modelR22.mat')
load('modelR22_FAD.mat')
load('modelR22_flavo.mat')
 
carbon_sources = {'EX_glc_D[e]'; 'EX_fru[e]'; 'EX_but[e]';
        'EX_hxa[e]'; 'EX_octa[e]'; 'EX_dca[e]'; 'EX_ddca[e]';
        'EX_ttdca[e]'; 'EX_hdca[e]'; 'EX_ocdca[e]'; 'EX_arach[e]';
        'EX_docosac[e]'; 'EX_lgnc[e]'; 'EX_hexc[e]';'EX_ala_L[e]';
        'EX_arg_L[e]'; 'EX_asn_L[e]'; 'EX_asp_L[e]'; 'EX_cys_L[e]';
        'EX_gln_L[e]'; 'EX_glu_L[e]'; 'EX_gly[e]'; 'EX_his_L[e]';
        'EX_ile_L[e]'; 'EX_leu_L[e]'; 'EX_lys_L[e]'; 'EX_met_L[e]';
        'EX_phe_L[e]'; 'EX_pro_L[e]'; 'EX_ser_L[e]'; 'EX_thr_L[e]';
        'EX_trp_L[e]'; 'EX_tyr_L[e]'; 'EX_val_L[e]'};
 
%set constraints
modelR22.lb(findRxnIDs(modelR22, carbon_sources))=-1;
modelR22_FAD.lb(findRxnIDs(modelR22_FAD, carbon_sources))=-1;
modelR22_flavo.lb(findRxnIDs(modelR22_flavo, carbon_sources))=-1;
modelR22.lb(findRxnIDs(modelR22, 'biomass_reaction'))=0.1;
modelR22_FAD.lb(findRxnIDs(modelR22_FAD, 'biomass_reaction'))=0.1;
modelR22_flavo.lb(findRxnIDs(modelR22_flavo, 'biomass_reaction'))=0.1;
%set objective function
modelR22 = changeObjective(modelR22, 'biomass_reaction');
modelR22_FAD = changeObjective(modelR22_FAD, 'biomass_reaction');
modelR22_flavo = changeObjective(modelR22_flavo, 'biomass_reaction');
optimizeCbModel(modelR22)
optimizeCbModel(modelR22_FAD)
optimizeCbModel(modelR22_flavo)
 
%introduce MADD knockout
gene = [];
gene{1} = '2110.1';
[modelR22_MADD,~,~,~] = deleteModelGenes(modelR22, gene);
[modelR22_FAD_MADD,~,~,~] = deleteModelGenes(modelR22_FAD, gene);
[modelR22_flavo_MADD,~,~,~] = deleteModelGenes(modelR22_flavo, gene);
%save temp_workspace
%reduce models
rModel_R = reduceModel(modelR22);
rModel_R_MADD = reduceModel(modelR22_MADD);
rModel_FAD = reduceModel(modelR22_FAD);
rModel_FAD_MADD = reduceModel(modelR22_FAD_MADD);
rModel_flavo = reduceModel(modelR22_flavo);
rModel_flavo_MADD = reduceModel(modelR22_flavo_MADD);
 
[rModel_R_MADD2,~,~,~] = deleteModelGenes(rModel_R, gene);
 
%save temp_workspace
%sample models using OptGpSampler
sModel_R = optGpSampler(rModel_R, [],10000, 2000, 4, 'glpk', 0);
sModel_R_MADD = optGpSampler(rModel_R_MADD2, [], 10000, 2000, 4, 'glpk', 0);
sModel_FAD = optGpSampler(rModel_FAD, [], 10000, 2000, 4, 'glpk', 0);
sModel_FAD_MADD = optGpSampler(rModel_FAD_MADD, [], 10000, 2000, 4, 'glpk', 0);
sModel_flavo = optGpSampler(rModel_flavo, [], 10000, 20000, 4, 'glpk', 0);
sModel_flavo_MADD = optGpSampler(rModel_flavo_MADD, [], 10000, 2000, 4, 'glpk', 0);
 
%% variable sampleMADD_R3D 
% columns 1-6 biomass_reaction, 6-12 ETF reaction, columns 12-18 total mFAO
% Column order in the set: 
% Recon3D, Recon3D_MADD, Recon3D_FAD, Recon3D_FAD_MADD, Recon3D_flavo, Recon3D_flaco_MADD 

rxnList = {'biomass_reaction', 'ETF'};
sampleMADD_R3D =[];
for i = 1:length(rxnList)
    id = findRxnIDs(sModel_R, rxnList(i));
    sampleMADD_R3D = [sampleMADD_R3D sModel_R.points(id, :)'];
    id = findRxnIDs(sModel_R_MADD, rxnList(i));
    if id == 0
        sampleMADD_R3D = [sampleMADD_R3D zeros(10000,1)];
    else
        sampleMADD_R3D = [sampleMADD_R3D sModel_R_MADD.points(id, :)'];
    end
    id = findRxnIDs(sModel_FAD, rxnList(i));
    sampleMADD_R3D = [sampleMADD_R3D sModel_FAD.points(id, :)'];
    id = findRxnIDs(sModel_FAD_MADD, rxnList(i));
    if id == 0
        sampleMADD_R3D = [sampleMADD_R3D zeros(10000,1)];
    else
        sampleMADD_R3D = [sampleMADD_R3D sModel_FAD_MADD.points(id, :)'];
    end
    id = findRxnIDs(sModel_flavo, rxnList(i));
    sampleMADD_R3D = [sampleMADD_R3D sModel_flavo.points(id, :)'];
    id = findRxnIDs(sModel_flavo_MADD, rxnList(i));
    if id == 0
        sampleMADD_R3D = [sampleMADD_R3D zeros(10000,1)];
    else
        sampleMADD_R3D = [sampleMADD_R3D sModel_flavo_MADD.points(id, :)'];
    end
end
 
% total mitochondrial FAO
AllETF_FAD = findRxnsFromMets(sModel_FAD, 'etfox[m]');
AllETF_flavo = findRxnsFromMets(sModel_flavo, 'etfox[m]');
Ctrl = [];
for i = 1:length(AllETF_FAD)
    id = findRxnIDs(sModel_R, AllETF_FAD(i));
    if id == 0
        Ctrl = [Ctrl; zeros(1,10000)];
    else
        Ctrl = [Ctrl; sModel_R.points(id, :)];
    end
end
 
Ctrl_MADD = [];
for i = 1:length(AllETF_FAD)
    id = findRxnIDs(sModel_R_MADD, AllETF_FAD(i));
    if id == 0
        Ctrl_MADD = [Ctrl_MADD; zeros(1,10000)];
    else
        Ctrl_MADD = [Ctrl_MADD; sModel_R_MADD.points(id, :)];
    end
end
 
FAD = [];
for i = 1:length(AllETF_FAD)
    id = findRxnIDs(sModel_FAD, AllETF_FAD(i));
    if id == 0
        FAD = [FAD; zeros(1,10000)];
    else
        FAD = [FAD; sModel_FAD.points(id, :)];
    end
end
 
FAD_MADD = [];
for i = 1:length(AllETF_FAD)
    id = findRxnIDs(sModel_FAD_MADD, AllETF_FAD(i));
    if id == 0
        FAD_MADD = [FAD_MADD; zeros(1,10000)];
    else
        FAD_MADD = [FAD_MADD; sModel_FAD_MADD.points(id, :)];
    end
end
 
 
FAD_flavo = [];
for i = 1:length(AllETF_flavo)
    id = findRxnIDs(sModel_flavo, AllETF_flavo(i));
    if id == 0
        FAD_flavo = [FAD_flavo; zeros(1,10000)];
    else
        FAD_flavo = [FAD_flavo; sModel_flavo.points(id, :)];
    end
end
FAD_flavo_MADD = [];
for i = 1:length(AllETF_flavo)
    id = findRxnIDs(sModel_flavo_MADD, AllETF_flavo(i));
    if id == 0
        FAD_flavo_MADD = [FAD_flavo_MADD; zeros(1,10000)];
    else
        FAD_flavo_MADD = [FAD_flavo_MADD; sModel_flavo_MADD.points(id, :)];
    end
end

sampleMADD_R3D = [sampleMADD_R3D sum(Ctrl,1)'];
sampleMADD_R3D = [sampleMADD_R3D sum(Ctrl_MADD,1)'];
sampleMADD_R3D = [sampleMADD_R3D sum(FAD,1)'];
sampleMADD_R3D = [sampleMADD_R3D sum(FAD_MADD,1)'];
sampleMADD_R3D = [sampleMADD_R3D sum(FAD_flavo,1)'];
sampleMADD_R3D = [sampleMADD_R3D sum(FAD_flavo_MADD,1)'];

