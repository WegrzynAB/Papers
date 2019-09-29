function geneIDs = findGenesFromEntrezIDs(model, genelist)

Genes = regexp(model.genes,'[\., or , (,)]','split');
modelGenes = cell(length(Genes),1);
for i = 1: length(Genes)
    modelGenes{i,1} = Genes{i,1}{1,1};
end

k = 1;
for i = 1:length(genelist) 
    ID = find(ismember(modelGenes,genelist{i}));
         if ~isempty(ID)
             if length(ID) > 1 
                 for j=1:length(ID)
                     geneIDs(k,1)  = model.genes(ID(j,1));
                     k = k+1;
                 end
             else
                geneIDs(k,1)  = model.genes(ID);
                k = k+1;
             end 
         end
end 
geneIDs(all(cellfun('isempty',geneIDs),2),:) = [];    
end
