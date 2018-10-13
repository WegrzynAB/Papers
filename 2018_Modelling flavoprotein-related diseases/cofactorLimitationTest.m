%% testing response of the system to the cofactorLimitation
File = 'modelR22_FAD.mat'; %'modelR3D_FAD.mat' for Recon3D 'modelR22_FAD' for Recon2.2
load(File);
%% Read Genes encoding enzymes that use  cofactor (provide correct pathway to your local S1_Table)
if strcmp('modelR3D_FAD.mat',File)
        [~, ~, raw] = xlsread('S1_Table.xlsx','Flavoproteins','G2:G112');
        ctrl = modelR3D_FAD;
        model = modelR3D_flavo;
    else
        [~, ~, raw] = xlsread('S1_Table.xlsx','Flavoproteins','J2:J112');
        ctrl = modelR22_FAD;
        model = modelR22_flavo;
end

% identify all flavoprotein-releted reactions in the models
[~,~,flavoProt_ctrl,~]=deleteModelGenes(ctrl,flavoprotMapped);
len_Ctrl = length(flavoProt_ctrl);
[~,~,flavoProt,~]=deleteModelGenes(model,flavoprotMapped);
len_Model = length(flavoProt);

% change objective to optimize for the flux of the flavoprotein-related
% reactions
model = changeObjective(model, flavoProt);
ctrl = changeObjective(ctrl, flavoProt_ctrl);
solution_Ctrl = optimizeCbModel(ctrl);
solution_Model = optimizeCbModel(model);

N=40;
gr = zeros(N+1,1);

ctrl = changeObjective(ctrl, 'FMNAT');
model = changeObjective(model, 'FMNAT');
sol_ctrl = optimizeCbModel(ctrl);
sol_model = optimizeCbModel(model);
maxCofactorFlux_ctrl = sol_ctrl.f;
maxCofactorFlux_model = sol_model.f;

model = changeObjective(model, flavoProt);
ctrl = changeObjective(ctrl, flavoProt_ctrl);

cofactorFlux = zeros(N+1,1);
cofactorFlux(1,1) = 0.02;

temp1 = ctrl;
temp2 = model;

Pos = 1;
for i = 0:0.2/N:0.2
    temp1=changeRxnBounds(temp1,'FMNAT',i,'b');
    temp2=changeRxnBounds(temp2,'FMNAT',i,'b');
    solution1 = optimizeCbModel(temp1,'max');
    solution2 = optimizeCbModel(temp2,'max');
    gr(Pos,1) = solution1.f/len_Ctrl;
    gr(Pos,2) = solution2.f/len_Model;
    cofactorFlux(Pos,1) = i;
    Pos = Pos+1;
end

clearvars -except modelR22_FAD modelR22_flavo gr cofactorFlux
