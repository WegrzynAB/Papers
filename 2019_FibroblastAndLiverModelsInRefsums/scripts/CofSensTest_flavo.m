%% testing sensitivity of the system to cofactor stoichiometry
initCobraToolbox(false);
File1 = 'E:/Dropbox/Cofactor Manuscript/scripts/models/modelR3D_FAD.mat'; %'modelR3D_FAD.mat' for Recon3D 'modelR22_FAD' for Recon2.2
load(File1);
File2 = 'E:/Dropbox/Cofactor Manuscript/scripts/models/modelR3D_flavo.mat';
load(File2);
File3 = 'E:/Dropbox/Cofactor Manuscript/scripts/models/modelR3D.mat';
load(File3);
modelOrig = 'Recon3';
%% Read Genes encoding enzymes that use  cofactor (provide correct pathway to your local S1_Table)
if strcmp('Recon3',modelOrig)
        [~, ~, raw] = xlsread('E:/Dropbox/Cofactor Manuscript/Additional_files/S1_Table.xlsx','Flavoproteins','G2:G112');
        ctrl = modelR3D_FAD;
        recon = Recon3DModel;
        flavo = modelR3D_flavo;
        model = modelR3D_X_c;
    else
        [~, ~, raw] = xlsread('E:/Dropbox/Cofactor Manuscript/Additional_files/S1_Table.xlsx','Flavoproteins','J2:J112');
        ctrl = modelR22_FAD;
        model = modelR22_X_c;
end

% identify all flavoprotein-releted reactions in the models
[~,~,flavoProt_ctrl,~]=deleteModelGenes(ctrl,flavoprotMapped);
len_Ctrl = length(flavoProt_ctrl);
[~,~,flavoProt_recon,~]=deleteModelGenes(recon,flavoprotMapped);
len_Recon = length(flavoProt_recon);
[~,~,flavoProt_flavo,~]=deleteModelGenes(flavo,flavoprotMapped);
len_Flavo = length(flavoProt_flavo);
[~,~,flavoProt_model,~]=deleteModelGenes(model,flavoprotMapped);
len_Model = length(flavoProt_model);

% change objective to optimize for the flux of the flavoprotein-related
% reactions
N=40;
gr = zeros(N+1,1);

ctrl = changeObjective(ctrl, 'FMNAT');
recon = changeObjective(recon, 'FMNAT');
flavo = changeObjective(flavo, 'FMNAT');
model = changeObjective(model, 'FMNAT');
sol_ctrl = optimizeCbModel(ctrl);
sol_recon = optimizeCbModel(recon);
sol_flavo = optimizeCbModel(flavo);
sol_model = optimizeCbModel(model);
maxCofactorFlux_ctrl = sol_ctrl.f;
maxCofactorFlux_recon = sol_recon.f;
maxCofactorFlux_flavo = sol_flavo.f;
maxCofactorFlux_model = sol_model.f;

ctrl = changeObjective(ctrl, flavoProt_ctrl);
recon = changeObjective(recon, flavoProt_recon);
flavo = changeObjective(flavo, flavoProt_flavo);
model = changeObjective(model, flavoProt_model);

cofactorFlux = zeros(N*3,1);

temp1 = ctrl;
temp2 = recon;
temp3 =  flavo;
temp4 = model;
Pos = 1;
tic
for i = 0:1.5/(N*3):1.5
    temp1=changeRxnBounds(temp1,'FMNAT',i,'b');
    temp2=changeRxnBounds(temp2,'FMNAT',i,'b');
    temp3=changeRxnBounds(temp3,'FMNAT',i,'b');
    temp4=changeRxnBounds(temp4,'FMNAT',i,'b');
    solution1 = optimizeCbModel(temp1,'max');
    solution2 = optimizeCbModel(temp2,'max');
    solution3 = optimizeCbModel(temp3,'max');
    solution4 = optimizeCbModel(temp4,'max');
    gr(Pos,1) = solution1.f/len_Ctrl;
    gr(Pos,2) = solution2.f/len_Recon;
    gr(Pos,3) = solution3.f/len_Flavo;
    gr(Pos,4) = solution4.f/len_Model;
    cofactorFlux(Pos,1) = i;
    Pos = Pos+1;
end
toc
clearvars -except modelR3D_FAD modelR3D_X_c modelR22_FAD modelR22_X_c gr cofactorFlux flavoprotMapped