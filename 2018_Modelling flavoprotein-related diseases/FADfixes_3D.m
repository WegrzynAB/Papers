%load('Recon22.mat', 'ReconR22')
%modelR22 = Recon22;
model = removeUnusedGenes(modelR22);
%% GPR fix
[~, ~, raw] = xlsread('/Users/Dave/Dropbox/Cofactor Manuscript/Additional_files/S2_Table_reviewed.xlsx','GPRdel','A2:A32');
%[~, ~, raw] = xlsread('D:\dropbox\Dropbox\Cofactor Manuscript\Additional_files\S2_Table_reviewed.xlsx','GPRdel','A2:A32');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,1);
Delete = cellVectors(:,1);
Del_3D = Delete(ismember(Delete,Recon3DModel.rxns));
%model = removeRxns(model, Delete);
model = removeRxns(model, Del_3D);
%Import the data to change GPR
[~, ~, raw0_0] = xlsread('/Users/Dave/Dropbox/Cofactor Manuscript/Additional_files/S2_Table_reviewed.xlsx','GPRfixes','A2:A289');
[~, ~, raw0_1] = xlsread('/Users/Dave/Dropbox/Cofactor Manuscript/Additional_files/S2_Table_reviewed.xlsx','GPRfixes','G2:H289');
[~, ~, raw0_2] = xlsread('/Users/Dave/Dropbox/Cofactor Manuscript/Additional_files/S2_Table_reviewed.xlsx','GPRfixes','I2:I289');
raw = [raw0_0,raw0_1,raw0_2];
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[1,4]);
raw = raw(:,[2,3]);
data = reshape([raw{:}],size(raw));

RxnID = cellVectors(:,1);
RxnID_3D = RxnID(ismember(RxnID, model.rxns));
LB = data(:,1);
UB = data(:,2);
GPRnew = cellVectors(:,2);

for i=1:length(RxnID)
    if ismember(RxnID(i), model.rxns)
        model = changeGeneAssociation(model,RxnID(i),GPRnew{i});
    end
end
model = changeRxnBounds(model,RxnID,LB,'l');
model = changeRxnBounds(model,RxnID,UB,'u');
modelR3D_GPR = model;
%% FAS fix 
% Import the data
[~, ~, raw] = xlsread('/Users/Dave/Dropbox/Cofactor Manuscript/Additional_files/S2_Table_reviewed.xlsx','FAS', 'A2:F15');
%[~, ~, raw] = xlsread('D:\dropbox\Dropbox\Cofactor Manuscript\Additional_files\S2_Table_reviewed.xlsx','FAS curation','A2:F15');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[1,2,3,4]);
raw = raw(:,[5,6]);
data = reshape([raw{:}],size(raw));

%Allocate imported array to column variable names
RxnFas = cellVectors(:,1);
FormulaFas = cellVectors(:,2);
GPRfas = cellVectors(:,3);
SubSysFas = cellVectors(:,4);
LBfas = data(:,1);
UBfas = data(:,2);

%model = modelR22_GPR;

for i=1:length(RxnFas)
        model = addReaction(model,RxnFas{i},FormulaFas{i},'','',LBfas(i),UBfas(i),0, SubSysFas{i},GPRfas{i});
end

modelR3D_FAS = model;
%% Perox fix
[~, ~, raw] = xlsread('/Users/Dave/Dropbox/Cofactor Manuscript/Additional_files/S2_Table_reviewed.xlsx','PeroxDel','B2:B20');
%[~, ~, raw] = xlsread('D:\dropbox\Dropbox\Cofactor Manuscript\Additional_files\S2_Table_reviewed.xlsx','PeroxDel','B2:B20');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,1);

Delete = cellVectors(:,1);

model = modelR3D_FAS;
model = removeRxns(model, Delete);

%Import the data to curate peroxisomal reactions
[~, ~, raw] = xlsread('/Users/Dave/Dropbox/Cofactor Manuscript/Additional_files/S2_Table_reviewed.xlsx','PeroxFAO');
raw = raw(2:end,[2,4,6,9:10]);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[1,2,3]);
raw = raw(:,[4,5]);
data = reshape([raw{:}],size(raw));

RxnID = cellVectors(:,1);
NewFormula = cellVectors(:,2);
NewGPR = cellVectors(:,3);
NewLB = data(:,1);
NewUB = data(:,2);


for i=1:length(RxnID)
    model = addReaction(model,RxnID{i},NewFormula{i},'','',NewLB(i),NewUB(i),0, 'Peroxisomal transport',NewGPR{i});
end

for i=1:length(RxnID)
   model = changeGeneAssociation(model,RxnID(i),NewGPR{i});
end
model = changeRxnBounds(model,RxnID,NewLB,'l');
model = changeRxnBounds(model,RxnID,NewUB,'u');
modelR3D_Per = model;
%% FAD fix
model = modelR3D_Per;
%Import the data
[~, ~, raw] = xlsread('/Users/Dave/Dropbox/Cofactor Manuscript/Additional_files/S2_Table_reviewed.xlsx','FAD delete','A2:A33');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,1);
delete = cellVectors(:,1);
%delete reactions and add FAD demand reaction
for i=1:length(delete)
    model = removeRxns(model, delete{i});
end
model = removeRxns(model, Del_3D);

%Import the data
[~, ~, raw] = xlsread('/Users/Dave/Dropbox/Cofactor Manuscript/Additional_files/S2_Table_reviewed.xlsx','FAD','A2:F161');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[1,2,3,6]);
raw = raw(:,[4,5]);
data = reshape([raw{:}],size(raw));

Abbreviation = cellVectors(:,1);
Formula = cellVectors(:,2);
Pathway = cellVectors(:,3);
LowerBound = data(:,1);
UpperBound = data(:,2);
GPR = cellVectors(:,4);
%add reactions
for i=1:length(Abbreviation)
    if ismember(Abbreviation(i),FadRxns_3D)
        model = addReaction(model,Abbreviation{i},Formula{i});
    end
end
for i=1:length(Fad_3Dfix)
        model = addReaction(model,Fad_3Dfix{i},Fad_3Dfix_form{i});
end

% for i=1:length(Abbreviation)
%     model = changeGeneAssociation(model,Abbreviation(i),GPR{i});
% end
% model = changeRxnBounds(model,Abbreviation,LowerBound,'l');
% model = changeRxnBounds(model,Abbreviation,UpperBound,'u');
% 
% model = addDemandReaction(model, 'fad[c]');
modelR3D_FAD = model;
%% Directionality fix
[~, ~, raw] = xlsread('/Users/Dave/Dropbox/Cofactor Manuscript/Additional_files/S2_Table_reviewed.xlsx','Directionality');
raw = raw(2:end,[1,5:6]);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,1);
raw = raw(:,[2,3]);
data = reshape([raw{:}],size(raw));
RxnID = cellVectors(:,1);
Lb = data(:,1);
Ub = data(:,2);

model = modelR3D_FAD;

model = changeRxnBounds(model,RxnID,Lb,'l');
model = changeRxnBounds(model,RxnID,Ub,'u');

model = removeUnusedGenes(model);
modelR3D_FAD = model;

%% Clear temporary variables
clearvars -except modelR204 modelR22 modelR22_FAD