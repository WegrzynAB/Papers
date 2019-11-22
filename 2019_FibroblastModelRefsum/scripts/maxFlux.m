function [flux, objective, carbon_source, normoxic] = maxFlux(model_org, carbon_source, objective, normoxic, media, verbose)

% Block import reactions:
model = blockAllImports(model_org);

% Define carbon source:
model = changeRxnBounds(model, carbon_source, -1, 'b');

% Define media:
model = changeRxnBounds(model, media, -1000, 'l');
model = changeRxnBounds(model, 'EX_co2[e]', 1000, 'u');
model = changeRxnBounds(model, 'EX_gthrd[e]', -1, 'l');
model = changeRxnBounds(model, 'EX_pnto_R[e]', -1, 'l');

if normoxic
    model = changeRxnBounds(model, 'EX_o2[e]', -1000, 'l');
end

% Specify objective and maximise:
model = changeObjective(model, objective);
FBAsolution = optimizeCbModel(model,'max','zero');
modelName = model.description;

if length(FBAsolution.x)>0
    flux = FBAsolution.f;
    %         result_filename_stem = strcat(modelName,'_',carbon_source);
    %         writeResult(model, FBAsolution, strcat(result_filename_stem, '.xls'));
    %         writeFusionTable(model, FBAsolution, strcat(result_filename_stem, '.txt'));
    if verbose == 1
        fprintf('%s\t%s\t%d\t%.2f\n', objective, carbon_source, normoxic, flux);
    end
else
    if verbose == 1
        fprintf('%s\t%s\t%d\tInfeasible\n', objective, carbon_source, normoxic);
    end
    flux = -1;
end

end