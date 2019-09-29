genes = {'6120.1';'22934.1';'6535.1';'686.1';'6888.1';...
    '6897.1';'7086.1';'875.1';'16.1';'79944.1';'10667.1';'57176.1';'18.1';...
    '55157.1';'1615.1';'728294.1';'1719.1';'124454.1';'7915.1';'2593.1';...
    '2628.1';'23438.1';'4524.1';'443.1';'501.1';'5917.1';'80704.1';'2108.1';...
    '2109.1';'2110.1';'35.1';'34.1';'33.1';'37.1';'3712.1';'2639.1';'36.1'};
listRxns = findRxnsFromGenes(fibroblast, genes);
fn = fieldnames(listRxns);
res = [];
for i=1:length(fn)
    temp = listRxns.(fn{i});
    for j = 1:length(temp(:,1))    
        rxn = listRxns.(fn{i}){j,1};
        sol_struct = optimizeCbModel(changeObjective(fibroblast, rxn));
        sol(j,1) = {fn{i}};
        sol(j,2) = {rxn};
        sol(j,3) = {sol_struct.f};
    end
    res = [res; sol];
    sol = cell(1,3);
end

