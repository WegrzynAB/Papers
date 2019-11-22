function [rxnList, geneIdRxns, IDs] = findRxnsFromEntrezIDs(model, genelist)

Genes = regexp(model.genes,'[\., or , (,)]','split');
modelGenes = cell(length(Genes),1);
for i = 1: length(Genes)
    modelGenes{i,1} = Genes{i,1}{1,1};
end

rxnList = {};
geneIdRxns = cell(length(genelist),4);
k = 1;
IDs = [];
for i = 1:length(genelist) 
    ID = find(ismember(modelGenes,genelist{i}));
         if ~isempty(ID)
             if ID > 1 
                 for j = 1:length(ID)
                    IDs = [IDs; ID(j,1)];
                 end
                 ID = ID(1,1);
             else
                 IDs = [IDs; ID];
             end 
             geneIdRxns(k,1)  = genelist(i);
             geneIdRxns(k,2)  = {ID};
             for j = 1:length(ID)
                 geneIdRxns(k,3) = {model.rxns(logical(model.rxnGeneMat(:,ID(j))))};
                 geneIdRxns(k,4) = {printRxnFormula(model,model.rxns(logical(model.rxnGeneMat(:,ID(j)))), false)}; 
                 rxnList = [rxnList; geneIdRxns{k,3}];
             end
             k = k+1;
         end
end 
geneIdRxns(all(cellfun('isempty',geneIdRxns),2),:) = [];
rxnList = unique(rxnList);      
end
