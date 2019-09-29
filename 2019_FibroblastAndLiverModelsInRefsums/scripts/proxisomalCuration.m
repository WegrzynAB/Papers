model = removeUnusedGenes(modelR3D_FAD);
%% Import the data for removed reactions
[~, ~, raw] = xlsread('E:\Dropbox\Refsum\manuscript\supplementary tables\ST1_peroxisomalFix.xlsx','Recon3D_del');
raw = raw(2:end,1);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,1);

% Allocate imported array to column variable names
RxnID_del = cellVectors(:,1);

% Clear temporary variables
clearvars raw cellVectors;

% delete reactions
model = removeRxns(model, RxnID_del);

%% Import the data for modified/added reactions
[~, ~, raw] = xlsread('E:\Dropbox\Refsum\manuscript\supplementary tables\ST1_peroxisomalFix.xlsx','Fixes');
raw = raw(2:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[1,2,5,6]);
raw = raw(:,[3,4]);
% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
RxnID = cellVectors(:,1);
Formula = cellVectors(:,2);
LB = data(:,1);
UB = data(:,2);
SubSystem = cellVectors(:,3);
GeneID = cellVectors(:,4);

% Clear temporary variables
clearvars data raw cellVectors;
%% Import metabolite information
opts = spreadsheetImportOptions("NumVariables", 4);
% Specify sheet and range
opts.Sheet = "Added_mets";
opts.DataRange = "A2:D23";
% Specify column names and types
opts.VariableNames = ["mets", "Name", "FormulaMet", "InChI"];
opts.SelectedVariableNames = ["mets", "Name", "FormulaMet", "InChI"];
opts.VariableTypes = ["string", "string", "string", "string"];
opts = setvaropts(opts, [1, 2, 3, 4], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2, 3, 4], "EmptyFieldRule", "auto");
% Import the data
tbl = readtable("E:\Dropbox\Refsum\manuscript\supplementary tables\ST1_peroxisomalFix.xlsx", opts, "UseExcel", false);
%Convert to output type
mets = tbl.mets;
Name = tbl.Name;
FormulaMet = tbl.FormulaMet;
InChI = tbl.InChI;
% Clear temporary variables
clear opts tbl
%%
% add/modify reactions

for i=1:length(RxnID)
        model = addReaction(model,RxnID{i},Formula{i},'','',LB(i),UB(i),0, SubSystem{i}, GeneID{i});
end
model = addExchangeRxn(model, 'but[e]');
model = addExchangeRxn(model, 'dca[e]');
model = checkDuplicateRxn(model,'S',1,0);
model.ub(findRxnIDs(model, 'EX_prist[e]')) = 0;
model.ub(findRxnIDs(model, 'EX_dmhptcrn[e]')) = 0;
model = checkDuplicateRxn(model, 'S', 1, 0);
model = removeUnusedGenes(model);

%update metabolite information
for i=1:length(mets)
    ID = findMetIDs(model, mets{i});
    model.metNames(ID) = {Name{i}};
    model.metFormulas(ID) = {FormulaMet{i}};
    model.metInChIString(ID) = {InChI{i}};
end
%% save final model with curated metabolism of phytanic acid
modelR3D_FAD_X = model;
save('modelR3D_FAD_X.mat', 'modelR3D_FAD_X')
