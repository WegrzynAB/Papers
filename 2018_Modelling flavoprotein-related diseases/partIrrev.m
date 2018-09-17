% Programmer:     Rienk Rienksma
% Last uptdated:  31 Aug 2012
% status:         Functional

% Function that can make a part of a model irreversible

% Input:
% model           a constraint-based model with a COBRA structure

% Optional input:
% rxnlist         a string of the names of the reactions that should be
%                 converted to irreversible reactions. A boolean vector,
%                 corresponding to the reactions in model.rxns, is
%                 also possible. This string or vector can also
%                 contain reactions that are already irreversible. Use the
%                 names in model.rxns (Default = all reactions in the
%                 model)

% Output
% irrevmodel      a model containing a number of irreversible reactions
%                 (depending on 'rxnlist'). The forward reactions have the
%                 extension '__fw', the backward reactions '__bw' and the
%                 reactions where both the bounds and the directionality
%                 are reversed '__rv'
% rxnlistnew      a boolean vector indicating all reactions in the new
%                 model

function [irrevmodel,rxnlistnew] = partIrrev(model,rxnlist)

if nargin < 2
    rxnlist = ones(length(model.rxns),1);
end

% Check is the input is a cell and convert to a boolean
if iscell(rxnlist)
    h = waitbar(0,'locating reactions in model...');
    rxnbool = zeros(length(model.rxns),1);
    for i = 1:length(model.rxns)
        waitbar(i/length(model.rxns))
        for j = 1:length(rxnlist)
            if strcmp(model.rxns(i),rxnlist(j))
                rxnbool(i,1) = 1;
            end
        end
    end
    if sum(rxnbool) ~= length(rxnlist)
        warning('Not all reactions are in the model, or there are duplicate reactions in the reaction list')
    end
    close(h)
    rxnbool = logical(rxnbool);
end

% Check is rxnlist is a vector
if isvector(rxnlist) && iscell(rxnlist) == 0
    rxnbool = logical(rxnlist);
end

% Check the reversibility of the reactions
irxns = zeros(length(model.rxns),1);
rrxns = zeros(length(model.rxns),1);
modelold = model;
x=1;

f = waitbar(0,'Creating irreversible model...');
for i = 1:length(rxnbool)
    waitbar(i/length(rxnbool))
    if rxnbool(i,1) == 1 && modelold.lb(i) == -1000 && modelold.ub(i) == 1000 % Identify the reversible reactions
        irxns(i+x-1,1)           = 1;
        irxns                    = [irxns(1:i+x-1,1);1;irxns(i+x:end,1)];
        rrxns                    = [rrxns(1:i+x-1,1);0;rrxns(i+x:end,1)];
        model.rxns               = [model.rxns(1:i+x-1,1);strcat(model.rxns(i+x-1,1),'__bw');model.rxns(i+x:end,1)]; % Add a backward reaction
        model.rxnNames           = [model.rxnNames(1:i+x-1,1);strcat(model.rxnNames(i+x-1,1),' backward');model.rxnNames(i+x:end,1)];
        model.rxns(i+x-1,1)      = strcat(model.rxns(i+x-1,1),'__fw'); % Change the reaction in a forward reaction
        model.rxnNames(i+x-1,1)  = strcat(model.rxnNames(i+x-1,1),' forward');
        model.S                  = [model.S(:,1:i+x-1) -model.S(:,i+x-1) model.S(:,i+x:end)]; %Change the Stoichiometric matrix the 'minus' sign changes the direction of the reaction
        model.lb(i+x-1,1)        = 0; % Set the lower bound of the forward reaction
        model.ub(i+x-1,1)        = 1000; % Set the upper bound of the forward reaction
        model.rev(i+x-1,1)       = 0; % Make the forward reaction irreversible
        model.rev                = [model.rev(1:i+x-1,1);0;model.rev(i+x:end,1)]; % Make the reverse reaction irreversible
        model.lb                 = [model.lb(1:i+x-1,1);0;model.lb(i+x:end,1)]; % Add the lower bound of the reverse reaction
        model.ub                 = [model.ub(1:i+x-1,1);1000;model.ub(i+x:end,1)]; % Add the upper bound of the reverse reaction
        model.c                  = [model.c(1:i+x-1,1);0;model.c(i+x:end,1)]; % Extend the c vector with a zero
        if isfield(model,'rules')
            model.rules = [model.rules(1:i+x-1,1);model.rules(i+x-1,1);model.rules(i+x:end,1)]; % Add rules to the reverse reaction
        end
        if isfield(model,'rxnGeneMat')
            model.rxnGeneMat = [model.rxnGeneMat(1:i+x-1,:);model.rxnGeneMat(i+x-1,:);model.rxnGeneMat(i+x:end,:)];
        end
        if isfield(model,'grRules')
            model.grRules = [model.grRules(1:i+x-1,1);model.grRules(i+x-1,1);model.grRules(i+x:end,1)];
        end
        if isfield(model,'subSystems')
            model.subSystems = [model.subSystems(1:i+x-1,1);model.subSystems(i+x-1,1);model.subSystems(i+x:end,1)];
        end
        if isfield(model,'confidenceScores')
            model.confidenceScores = [model.confidenceScores(1:i+x-1,1);model.confidenceScores(i+x-1,1);model.confidenceScores(i+x:end,1)];
        end
        if isfield(model,'rxnReferences')
            model.rxnReferences = [model.rxnReferences(1:i+x-1,1);model.rxnReferences(i+x-1,1);model.rxnReferences(i+x:end,1)];
        end
        if isfield(model,'rxnECNumbers')
            model.rxnECNumbers = [model.rxnECNumbers(1:i+x-1,1);model.rxnECNumbers(i+x-1,1);model.rxnECNumbers(i+x:end,1)];
        end   
        if isfield(model,'rxnNotes')
            model.rxnNotes = [model.rxnNotes(1:i+x-1,1);model.rxnNotes(i+x-1,1);model.rxnNotes(i+x:end,1)];
        end    
        if isfield(model,'rxnKeggID')
            model.rxnKeggID = [model.rxnKeggID(1:i+x-1,1);model.rxnKeggID(i+x-1,1);model.rxnKeggID(i+x:end,1)];
        end
        if isfield(model,'rxnConfidenceEcoIDA')
            model.rxnConfidenceEcoIDA = [model.rxnConfidenceEcoIDA(1:i+x-1,1);model.rxnConfidenceEcoIDA(i+x-1,1);model.rxnConfidenceEcoIDA(i+x:end,1)];
        end
        if isfield(model,'rxnConfidenceScores')
            model.rxnConfidenceScores = [model.rxnConfidenceScores(1:i+x-1,1);model.rxnConfidenceScores(i+x-1,1);model.rxnConfidenceScores(i+x:end,1)];
        end
        if isfield(model,'rxnsboTerm')
            model.rxnsboTerm = [model.rxnsboTerm(1:i+x-1,1);model.rxnsboTerm(i+x-1,1);model.rxnsboTerm(i+x:end,1)];
        end
        x = x+1; % Increase the loopcounter by 1
    end
    if rxnbool(i,1) == 1 && modelold.lb(i) == -1000 && modelold.ub(i) == 0 % Identify the reversed reactions
        rrxns(i+x-1,1)           = 1;
        model.rxns(i+x-1,1)      = strcat(model.rxns(i+x-1,1),'__rv'); % Change the reaction in a forward reaction
        model.rxnNames(i+x-1,1)  = strcat(model.rxnNames(i+x-1,1),' reversed');
        model.S                  = [model.S(:,1:i+x-2) -model.S(:,i+x-1) model.S(:,i+x:end)]; %Change the Stoichiometric matrix the 'minus' sign changes the direction of the reaction
        model.lb(i+x-1,1)        = 0; % Set the lower bound of the forward reaction
        model.ub(i+x-1,1)        = 1000; % Set the upper bound of the forward reaction
        model.rev(i+x-1,1)       = 0; % Make the forward reaction irreversible
    end
end
close(f)
irrevmodel = model;
rxnlistnew = irxns+rrxns;

g = waitbar(0,'Looking for unchanged reactions...');
rxnlist = modelold.rxns(rxnbool);
for i = 1:length(rxnlist)
    waitbar(i/length(rxnlist))
    for j = 1:length(model.rxns)
        if strcmp(rxnlist(i),model.rxns(j))
            rxnlistnew(j,1) = 1;
        end
    end
end
close(g)
rxnlistnew = logical(rxnlistnew);