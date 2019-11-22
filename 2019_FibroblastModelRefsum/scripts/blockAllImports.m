function modelClosed = blockAllImports(model)

modelexchanges1 = strmatch('Ex_',model.rxns);
modelexchanges4 = strmatch('EX_',model.rxns);
modelexchanges2 = strmatch('DM_',model.rxns);
modelexchanges3 = strmatch('sink_',model.rxns);
% also close biomass reactions
BM= (find(~cellfun(@isempty,strfind(lower(model.mets),'bioma'))));
selExc = (find( full((sum(abs(model.S)==1,1) ==1) & (sum(model.S~=0) == 1))))';

modelexchanges = unique([modelexchanges1;modelexchanges2;modelexchanges3;modelexchanges4;selExc;BM]);
model.lb(find(ismember(model.rxns,model.rxns(modelexchanges))))=0;
model.c = zeros(length(model.rxns),1);
model.ub(selExc)=1000;
modelClosed = model;

end