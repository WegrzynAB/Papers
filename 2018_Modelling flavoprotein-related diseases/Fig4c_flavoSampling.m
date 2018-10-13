load('modelR22_FAD')
load('modelR22_flavo')

ctrl = modelR22_FAD;
model = modelR22_flavo;
%% Read Genes encoding enzymes that use  cofactor
[~, ~, raw] = xlsread('S1_Table.xlsx','Flavoproteins','J2:J112');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,1);
flavoprot = cellVectors(:,1);
clearvars raw cellVectors;
%% sampling with the gradually reduced flux through FMNAT reaction
n=1;
for i = 1:1:length(flavoprot)
    for j = 1:1:length(model.genes)
        x = strcmp(flavoprot(i),model.genes(j));
        if x == 1
            flavoprotMapped(n,1) = model.genes(j);
            n=n+1;
        end
    end
end

[~,~,flavoProt_ctrl,~]=deleteModelGenes(ctrl,flavoprotMapped); 
[~,~,flavoProt,~]=deleteModelGenes(model,flavoprotMapped);

ctrl = changeObjective(ctrl, flavoProt_ctrl); 
model = changeObjective(model, flavoProt);
rModel_fad = reduceModel(ctrl);
rModel_flavo = reduceModel(model);
A = sum(ismember(flavoProt_ctrl, rModel_fad.rxns));
B = sum(ismember(flavoProt, rModel_flavo.rxns));
sumC = [];
sumF = [];

for i=0:0.12/20:0.12
    rModel_flavo = changeRxnBounds(rModel_flavo, 'FMNAT', i, 'b');
    rModel_fad = changeRxnBounds(rModel_fad, 'FMNAT', i, 'b');
    sModel_fad = optGpSampler(rModel_fad, [],10000, 2000, 4, 'glpk', 0);
    sModel_flavo = optGpSampler(rModel_flavo, [],10000, 2000, 4, 'glpk', 0);
    fadFlux = zeros(1,10000);
    for j=1:length(flavoProt_ctrl)
        ID = findRxnIDs(sModel_fad, flavoProt_ctrl(j));
        if ID == 0
            fadFlux(j,:) = zeros(1,10000);
        else
            fadFlux(j,:) = sModel_fad.points(ID,:);
        end
    end
    sumC = [sumC sum(fadFlux)'];
    
    for k=1:length(flavoProt_ctrl)
    ID = findRxnIDs(sModel_flavo, flavoProt_ctrl(k));
        if ID == 0
            fadFlux(k,:) = zeros(1,10000);
        else
            fadFlux(k,:) = sModel_flavo.points(ID,:);
        end
    end
    sumF = [sumF sum(fadFlux)'];
end

avrC = mean(sumC./A)';
stdC = std(sumC./A)';
avrF = mean(sumF./B)';
stdF = std(sumF./B)';

flavoSampled = [avrC stdC avrF stdF];