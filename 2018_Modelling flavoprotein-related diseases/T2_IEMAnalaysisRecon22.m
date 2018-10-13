% Edited to work with Recon2.2 version by Agnieszka Wegrzyn, 2017
 
% This script requires
%  - the COBRA toolbox - https://github.com/opencobra/cobratoolbox
%  - IBM clpex or GLPK LP solver
% - fastFVA (Gudmunsson, S., Thiele, I., "Computationally efficient flux
% variability analysis", BMC Bioinf, 11:489 (2010).)
%
%
% - Please note that this script performs more than 400 fastFVA
% calculations, each of which takes at least 35 sec depending on computer
% configuration and the number of parallel nodes (default 4).
%
%
% Ines Thiele, http://thielelab.eu, 2013, 2017

initCobraToolbox;
%% define variables:
% define solver for fastFVA
solver ='glpk';
% give number of workers for parallelization of fastFVA
nworkers=4;
% define input file name
File ='Recon22model_flavo';

load(File);
modelRecon22 = modelR22_flavo;
modelRecon22Model = modelR22_flavo;

% set all uptakes to -1 and all secretions to 1000
modelRecon22Model.lb(strmatch('EX_',modelRecon22Model.rxns))=-1;
modelRecon22Model.ub(strmatch('EX_',modelRecon22Model.rxns))=1000;
clear R22;
R22.model = modelRecon22Model;

modelRecon22=findSExRxnInd(modelRecon22); %modelRecon22.ExchRxnBool is made
modelRecon22Model=findSExRxnInd(modelRecon22Model); %modelRecon22Model.ExchRxnBool
modelRecon22.EXRxnBool = modelRecon22.ExchRxnBool; %old version of function so these have to be reversed
modelRecon22Model.EXRxnBool = modelRecon22Model.ExchRxnBool; %old version of function so these have to be reversed
R22.model=findSExRxnInd(R22.model);
R22.model.EXRxnBool = R22.model.ExchRxnBool; %old version of function so these have to be reversed
%% Map IEMs onto genes
load('IEM_compendiumR22');
% unique genes
R22.UniqueGenes=unique(modelRecon22.genes);
R22.IEMGenes=IEMsR22(ismember(IEMsR22(:,2),R22.UniqueGenes),2);
R22.IEMNames=IEMsR22(ismember(IEMsR22(:,2),R22.UniqueGenes),1);
%% perform fastFVA calculations
SetWorkerCount(nworkers);

ExR22 = modelRecon22Model.rxns(modelRecon22Model.EXRxnBool);
ExR22ID = find(modelRecon22Model.EXRxnBool);
% find all genes that are in model
[R22.NonuniqueGenes,rem]=strtok(modelRecon22Model.genes,'.');
% delete IEM gene
for i =1 :length(R22.IEMGenes)
    i
    tmp = strmatch(R22.IEMGenes(i),R22.NonuniqueGenes,'exact');
    Genes = modelRecon22Model.genes(tmp);
    [modelR22IEM,hasEffect,constrRxnNames,deletedGenes] = deleteModelGenes(modelRecon22Model,Genes);
    if hasEffect
        % test IEM case
            tic;[R22.IEMs.(char(strrep((strcat('G_',R22.IEMGenes(i))),':', '_'))).Disease.minFlux,R22.IEMs.(char(strrep((strcat('G_',R22.IEMGenes(i))),':', '_'))).Disease.maxFlux] = fastFVA(modelR22IEM,0,'max', solver);toc
        % test healthy case for IEM gene - healthy state --> force flux through
        % all gene asso rxns
        [modelR22Up,hasEffect,constrRxnNames,upregulatedGenes] = upRegulateModelGenes(modelRecon22Model,Genes,0.05);
        [R22.IEMs.(char(strrep((strcat('G_',R22.IEMGenes(i))),':', '_'))).Healthy.minFlux,R22.IEMs.(char(strrep((strcat('G_',R22.IEMGenes(i))),':', '_'))).Healthy.maxFlux] = fastFVA(modelR22Up,0,'max', solver);
    else
        R22.IEMEffectTestSolution(i,1)=-1;
    end
    save R22Results_tmp0 R22
end

%% % data analysis IEM
tol = 1e-6;
cnt= 1;
Factor= 0.99;

clear UU DU DD UD cnt i j R*BDone  R*Biomarker R*Change Change Biomarker R*ShlomiOmim WW tmp c ans NBB ExR* FR* Omim R*Genes

%%R2
clear R22Biomarker
FR22= fieldnames(R22.IEMs);
ExR22 = find(R22.model.EXRxnBool);

%ExR2=find(ismember(R2.model.rxns,R2.IEMs.(FR2{1}).Biomarkers));
cnt= 1;
for i = 1 : length(FR22)
    R22.IEMs.(FR22{i}).Biomarker=[];
    R22.IEMs.(FR22{i}).Healthy.minFlux(abs(R22.IEMs.(FR22{i}).Healthy.minFlux)<=tol)=0;
    R22.IEMs.(FR22{i}).Healthy.maxFlux(abs(R22.IEMs.(FR22{i}).Healthy.maxFlux)<=tol)=0;
    R22.IEMs.(FR22{i}).Disease.minFlux(abs(R22.IEMs.(FR22{i}).Disease.minFlux)<=tol)=0;
    R22.IEMs.(FR22{i}).Disease.maxFlux(abs(R22.IEMs.(FR22{i}).Disease.maxFlux)<=tol)=0;
    for j = 1 : length(ExR22)
        a1 = R22.IEMs.(FR22{i}).Healthy.minFlux(ExR22(j));
        a2 = R22.IEMs.(FR22{i}).Healthy.maxFlux(ExR22(j));
        b1 = R22.IEMs.(FR22{i}).Disease.minFlux(ExR22(j));
        b2 = R22.IEMs.(FR22{i}).Disease.maxFlux(ExR22(j));

        % disease secreted
        if a2 < Factor*b1
            % AAAAAAAAAA
            %            BBBBBBBBBB
            R22.IEMs.(FR22{i}).Biomarker(j,1) = 2;%strong
            R22.IEMsList(cnt,1)= FR22(i);
            R22.IEMsList(cnt,2)= R22.model.rxns(ExR22(j));
            R22.IEMsList(cnt,3)= num2cell(R22.IEMs.(FR22{i}).Biomarker(j,1));
            R22.IEMsList(cnt,4)= num2cell(a1);
            R22.IEMsList(cnt,5)= num2cell(a2);
            R22.IEMsList(cnt,6)= num2cell(b1);
            R22.IEMsList(cnt,7)= num2cell(b2);
            cnt = cnt +1;
            % disease more taken up
        elseif b2 < Factor*a1
            %            AAAAAAAAAA
            % BBBBBBBBBB
            R22.IEMs.(FR22{i}).Biomarker(j,1) = -2;%strong
            R22.IEMsList(cnt,1)= FR22(i);
            R22.IEMsList(cnt,2)= R22.model.rxns(ExR22(j));
            R22.IEMsList(cnt,3)= num2cell(R22.IEMs.(FR22{i}).Biomarker(j,1));
            R22.IEMsList(cnt,4)= num2cell(a1);
            R22.IEMsList(cnt,5)= num2cell(a2);
            R22.IEMsList(cnt,6)= num2cell(b1);
            R22.IEMsList(cnt,7)= num2cell(b2);
            cnt = cnt +1;

        elseif b1>=0 && a1 <= b1 && a2 < Factor*b2
            %      AAAAAAAAAA
            %       0  BBBBBBBBBB
            R22.IEMs.(FR22{i}).Biomarker(j,1) = 1;%weak
            R22.IEMsList(cnt,1)= FR22(i);
            R22.IEMsList(cnt,2)= R22.model.rxns(ExR22(j));
            R22.IEMsList(cnt,3)=num2cell(R22.IEMs.(FR22{i}).Biomarker(j,1));
            R22.IEMsList(cnt,4)= num2cell(a1);
            R22.IEMsList(cnt,5)= num2cell(a2);
            R22.IEMsList(cnt,6)= num2cell(b1);
            R22.IEMsList(cnt,7)= num2cell(b2);
            cnt = cnt +1;
        elseif b1>=0 && a1 < Factor*b1 && a2 <= b2
            %       AAAAAAAAAA
            %       0  BBBBBBBBBB
            R22.IEMs.(FR22{i}).Biomarker(j,1) = 1;%weak
            R22.IEMsList(cnt,1)= FR22(i);
            R22.IEMsList(cnt,2)= R22.model.rxns(ExR22(j));
            R22.IEMsList(cnt,3)= num2cell(R22.IEMs.(FR22{i}).Biomarker(j,1));
            R22.IEMsList(cnt,4)= num2cell(a1);
            R22.IEMsList(cnt,5)= num2cell(a2);
            R22.IEMsList(cnt,6)= num2cell(b1);
            R22.IEMsList(cnt,7)= num2cell(b2);
            cnt = cnt +1;
        elseif a1>=0 && b1 < Factor*a1 && b2 <= a2
            %  0    AAAAAAAAAA
            %     BBBBBBBBB
            R22.IEMs.(FR22{i}).Biomarker(j,1) = -1;%weak
            R22.IEMsList(cnt,1)= FR22(i);
            R22.IEMsList(cnt,2)= R22.model.rxns(ExR22(j));
            R22.IEMsList(cnt,3)= num2cell(R22.IEMs.(FR22{i}).Biomarker(j,1));
            R22.IEMsList(cnt,4)= num2cell(a1);
            R22.IEMsList(cnt,5)= num2cell(a2);
            R22.IEMsList(cnt,6)= num2cell(b1);
            R22.IEMsList(cnt,7)= num2cell(b2);
            cnt = cnt +1;
        elseif a1>=0 && b1 <= a1 && b2 < Factor*a2
            %  0  AAAAAAAAAA
            % BBBBBBBBB
            R22.IEMs.(FR22{i}).Biomarker(j,1) = -1;%weak
            R22.IEMsList(cnt,1)= FR22(i);
            R22.IEMsList(cnt,2)= R22.model.rxns(ExR22(j));
            R22.IEMsList(cnt,3)= num2cell(R22.IEMs.(FR22{i}).Biomarker(j,1));
            R22.IEMsList(cnt,4)= num2cell(a1);
            R22.IEMsList(cnt,5)= num2cell(a2);
            R22.IEMsList(cnt,6)= num2cell(b1);
            R22.IEMsList(cnt,7)= num2cell(b2);
            cnt = cnt +1;

        elseif b1<=0 && b2 >= 0 && Factor*a1 < b1 && a2 <= b2
            % AAAAAAAAAA
            %   BBBBB0BBBB
            R22.IEMs.(FR22{i}).Biomarker(j,1) = 1;%weak
            R22.IEMsList(cnt,1)= FR22(i);
            R22.IEMsList(cnt,2)= R22.model.rxns(ExR22(j));
            R22.IEMsList(cnt,3)= num2cell(R22.IEMs.(FR22{i}).Biomarker(j,1));
            R22.IEMsList(cnt,4)= num2cell(a1);
            R22.IEMsList(cnt,5)= num2cell(a2);
            R22.IEMsList(cnt,6)= num2cell(b1);
            R22.IEMsList(cnt,7)= num2cell(b2);
            cnt = cnt +1;
        elseif b1<=0 && b2 >= 0 && a1 <= b1 && a2 < Factor*b2
            % AAAAAAA0AAA
            %   BBBBB0BBBB
            R22.IEMs.(FR22{i}).Biomarker(j,1) = 1;%weak
            R22.IEMsList(cnt,1)= FR22(i);
            R22.IEMsList(cnt,2)= R22.model.rxns(ExR22(j));
            R22.IEMsList(cnt,3)= num2cell(R22.IEMs.(FR22{i}).Biomarker(j,1));
            R22.IEMsList(cnt,4)= num2cell(a1);
            R22.IEMsList(cnt,5)= num2cell(a2);
            R22.IEMsList(cnt,6)= num2cell(b1);
            R22.IEMsList(cnt,7)= num2cell(b2);
            cnt = cnt +1;

        elseif a1<=0 && a2>=0 && b1 <= a1 && b2 < Factor*a2
            %    AAAAA0AAAA
            %  BBBBBBB0BB
            R22.IEMs.(FR22{i}).Biomarker(j,1) = -1;%weak
            R22.IEMsList(cnt,1)= FR22(i);
            R22.IEMsList(cnt,2)= R22.model.rxns(ExR22(j));
            R22.IEMsList(cnt,3)=num2cell(R22.IEMs.(FR22{i}).Biomarker(j,1));
            R22.IEMsList(cnt,4)= num2cell(a1);
            R22.IEMsList(cnt,5)= num2cell(a2);
            R22.IEMsList(cnt,6)= num2cell(b1);
            R22.IEMsList(cnt,7)= num2cell(b2);
            cnt = cnt +1;
        elseif a1<=0 && a2>=0 && Factor*b1 < a1 && b2 <= a2
            %    AAAAA0AAAA
            %  BBBBBBB0BB
            R22.IEMs.(FR22{i}).Biomarker(j,1) = -1;%weak
            R22.IEMsList(cnt,1)= FR22(i);
            R22.IEMsList(cnt,2)= R22.model.rxns(ExR22(j));
            R22.IEMsList(cnt,3)=num2cell(R22.IEMs.(FR22{i}).Biomarker(j,1));
            R22.IEMsList(cnt,4)= num2cell(a1);
            R22.IEMsList(cnt,5)= num2cell(a2);
            R22.IEMsList(cnt,6)= num2cell(b1);
            R22.IEMsList(cnt,7)= num2cell(b2);
            cnt = cnt +1;

        elseif b2 <= 0 && Factor*a1 < b1 && a2 <= b2
            % AAAAAAAAAA   0
            %    BBBBBBBBB 0
            R22.IEMs.(FR22{i}).Biomarker(j,1) = 1;%weak
            R22.IEMsList(cnt,1)= FR22(i);
            R22.IEMsList(cnt,2)= R22.model.rxns(ExR22(j));
            R22.IEMsList(cnt,3)= num2cell(R22.IEMs.(FR22{i}).Biomarker(j,1));
            R22.IEMsList(cnt,4)= num2cell(a1);
            R22.IEMsList(cnt,5)= num2cell(a2);
            R22.IEMsList(cnt,6)= num2cell(b1);
            R22.IEMsList(cnt,7)= num2cell(b2);
            cnt = cnt +1;
        elseif a2<=0 && Factor*b1 < a1 && b2 <= a2
            %    AAAAAAAAA 0
            %  BBBBBBBBB   0
            R22.IEMs.(FR22{i}).Biomarker(j,1) = -1;%weak
            R22.IEMsList(cnt,1)= FR22(i);
            R22.IEMsList(cnt,2)= R22.model.rxns(ExR22(j));
            R22.IEMsList(cnt,3)=num2cell(R22.IEMs.(FR22{i}).Biomarker(j,1));
            R22.IEMsList(cnt,4)= num2cell(a1);
            R22.IEMsList(cnt,5)= num2cell(a2);
            R22.IEMsList(cnt,6)= num2cell(b1);
            R22.IEMsList(cnt,7)= num2cell(b2);
            cnt = cnt +1;
        elseif b2 <= 0 && a1 <= b1 && Factor*a2 < b2
            % AAAAAAAAAA   0
            %    BBBBBBBBB 0
            R22.IEMs.(FR22{i}).Biomarker(j,1) = 1;%weak
            R22.IEMsList(cnt,1)= FR22(i);
            R22.IEMsList(cnt,2)= R22.model.rxns(ExR22(j));
            R22.IEMsList(cnt,3)= num2cell(R22.IEMs.(FR22{i}).Biomarker(j,1));
            R22.IEMsList(cnt,4)= num2cell(a1);
            R22.IEMsList(cnt,5)= num2cell(a2);
            R22.IEMsList(cnt,6)= num2cell(b1);
            R22.IEMsList(cnt,7)= num2cell(b2);
            cnt = cnt +1;
        elseif a2<=0 && b1 <= a1 && Factor*b2 < a2
            %    AAAAAAAAA 0
            %  BBBBBBBBB   0
            R22.IEMs.(FR22{i}).Biomarker(j,1) = -1;%weak
            R22.IEMsList(cnt,1)= FR22(i);
            R22.IEMsList(cnt,2)= R22.model.rxns(ExR22(j));
            R22.IEMsList(cnt,3)=num2cell(R22.IEMs.(FR22{i}).Biomarker(j,1));
            R22.IEMsList(cnt,4)= num2cell(a1);
            R22.IEMsList(cnt,5)= num2cell(a2);
            R22.IEMsList(cnt,6)= num2cell(b1);
            R22.IEMsList(cnt,7)= num2cell(b2);
            cnt = cnt +1;

        else
            R22.IEMs.(FR22{i}).Biomarker(j,1) = 0;%none
            R22.IEMsList(cnt,1)= FR22(i);
            R22.IEMsList(cnt,2)= R22.model.rxns(ExR22(j));
            R22.IEMsList(cnt,3)= num2cell(R22.IEMs.(FR22{i}).Biomarker(j,1));
            R22.IEMsList(cnt,4)= num2cell(a1);
            R22.IEMsList(cnt,5)= num2cell(a2);
            R22.IEMsList(cnt,6)= num2cell(b1);
            R22.IEMsList(cnt,7)= num2cell(b2);
            cnt = cnt +1;
        end
    end
    R22Biomarker(:,i)=R22.IEMs.(FR22{i}).Biomarker;
end

R22.Biomarker=R22Biomarker;
% compounds that are no biomarker
a = 1 ;
for i = 1 :size(R22.Biomarker,1)
    if nnz(R22.Biomarker(i,:)) > 0
        R22.BiomarkerReal(a,:) = R22.Biomarker(i,:);
        R22.BiomarkerRealName(a,1) = R22.model.rxns(ExR22(i));
        a = a+1;
    end
end
% disease without biomarker

a = 1 ;
b=1;
R22.BiomarkerRealDisease=[];
for i = 1 :size(R22.Biomarker,2)
    if nnz(R22.Biomarker(:,i)) > 0
        R22.BiomarkerRealDisease(:,a) = R22.BiomarkerReal(:,i);
        R22.BiomarkerRealGenes(a,1) = regexprep(FR22(i),'G_','');
        a = a+1;
    else
        R22.NoBiomarkerGenes(b,1) = FR22(i);
        b = b +1;
    end
end

%IEMs without biomarker
R22.IEMsWOBiomarker= setdiff(R22.IEMGenes,R22.BiomarkerRealGenes);

R22Genes=regexprep(R22.IEMsList(:,1),'G_','');
% map swagatikas biomarkers and see what's happening
load IEM_biomarkerListR22;
if strcmp('2017_03_22_Recon3d_consistencyCheck',File)
    Swagatika_IEM_biomarker(:,2) = regexprep(Swagatika_IEM_biomarker(:,2),'\(','\[');
    Swagatika_IEM_biomarker(:,2) = regexprep(Swagatika_IEM_biomarker(:,2),'\)','\]');
    Swagatika_IEM_biomarker(:,2) = regexprep(Swagatika_IEM_biomarker(:,2),'EX_glc\[e\]','EX_glc_D[e]');
    Swagatika_IEM_biomarker(:,2) = regexprep(Swagatika_IEM_biomarker(:,2),'\(','\[');
    % check that all exchanges are matching
    Swagatika_IEM_biomarker(:,2) = strcat(Swagatika_IEM_biomarker(:,2),'[e]');
    Swagatika_IEM_biomarker(:,2) = regexprep(Swagatika_IEM_biomarker(:,2),'\[e\]\[e\]','\[e\]');
    setdiff(Swagatika_IEM_biomarker(:,2), R22.model.rxns(ExR22))
end

UU=0;DD=0;UD=0;DU=0;OO=0; OU=0; OD=0; NBB=0;
cnt = 1;
R22BDone = {};
R22BD={};   c = 1;
for i = 1 : size(Swagatika_IEM_biomarker,1)
    tmp = strmatch(Swagatika_IEM_biomarker(i,1),R22.IEMNames,'exact');

    if ~isempty(tmp)
        List=strmatch(R22.IEMGenes(tmp(1)),strrep(R22Genes,'_', ':'),'exact');
        Biomarker = Swagatika_IEM_biomarker(i,2);
        Change = str2num(Swagatika_IEM_biomarker{i,3});
        for j = 1 : length(List)
            R22Biomarker = R22.IEMsList(List(j),2);
            R22BD=strcat(R22Biomarker,'_',Swagatika_IEM_biomarker(i,1));
            R22Change = R22.IEMsList(List(j),3);

            if strcmp(Biomarker, R22Biomarker) && isempty(strmatch(R22BD,R22BDone,'exact'))
                % get change
                R22Change = R22.IEMsList{List(j),3};
                if R22Change ==2
                    R22Change =1;
                elseif R22Change== -2
                    R22Change=-1;
                end
                if Change ==R22Change
                    if Change==-1
                        DD = DD + 1;
                    elseif  Change==1
                        UU = UU +1;
                    elseif Change==0
                        OO = OO +1;
                    end
                elseif  Change==-1 && R22Change == 1
                    DU = DU + 1;
                elseif  Change==1 && R22Change == -1
                    UD = UD + 1;
                elseif Change == 0 && R22Change == 1
                    OU = OU +1;
                elseif Change == 0 && R22Change == -1
                    OD = OD +1;
                elseif R22Change == 0
                    NBB = NBB +1;
                end
                
                R22.Swagatika_IEM_biomarker(cnt,1) = Swagatika_IEM_biomarker(i);
                R22.Swagatika_IEM_biomarker(cnt,2) = Biomarker;
                R22.Swagatika_IEM_biomarker(cnt,3) = Swagatika_IEM_biomarker(i,3);
                R22.Swagatika_IEM_biomarker(cnt,4) = R22Biomarker;
                R22.Swagatika_IEM_biomarker(cnt,5) = R22.IEMsList(List(j),3);
                R22.Swagatika_IEM_biomarker(cnt,6) = R22.IEMsList(List(j),4);
                R22.Swagatika_IEM_biomarker(cnt,7) = R22.IEMsList(List(j),5);
                R22.Swagatika_IEM_biomarker(cnt,8) = R22.IEMsList(List(j),6);
                R22.Swagatika_IEM_biomarker(cnt,9) = R22.IEMsList(List(j),7);
                cnt = cnt +1;
                R22BDone(c)=R22BD;
                c=c+1;
            end
        end
    end
end

R22.SwagatikaIEM_biomarker.UU =UU;
R22.SwagatikaIEM_biomarker.DD =DD;
R22.SwagatikaIEM_biomarker.UD =UD;
R22.SwagatikaIEM_biomarker.DU =DU;
R22.SwagatikaIEM_biomarker.OO =OO;
R22.SwagatikaIEM_biomarker.OU =OU;
R22.SwagatikaIEM_biomarker.OD =OD;
R22.SwagatikaIEM_biomarker.NBB =NBB;
R22.SwagatikaIEM_biomarker.Sensitivity = UU/(UU+UD); %Recall
R22.SwagatikaIEM_biomarker.Specificity = DD/(DD+DU); %True negative rate
R22.SwagatikaIEM_biomarker.Precision = UU/(UU+DU);
R22.SwagatikaIEM_biomarker.Accuracy = (DD+UU)/(UU+DD+UD+DU);
%[p,x2] = chisquarecont([UU UD; DU DD]);
%R2.SwagatikaIEM_biomarker.chiSq = ((UU*DD-DU*UD)^2*(UU+DD+UD+DU))/((UU+UD)*(DU+DD)*(UU+DU)*(UD+DD));
% calculates the hypergeometric p value
%R2.SwagatikaIEM_biomarker.FischerExact= (factorial(vpi(UU + UD))* factorial(vpi(DD + DU))* factorial(vpi(UU + DU))* factorial(vpi(DD + UD)))/(factorial(vpi(UU))*factorial(vpi(DD))*factorial(vpi(UD))*factorial(vpi(DU))*factorial(vpi(UU+DD+UD+DU)));
R22.SwagatikaIEM_biomarker.FischerExact= (((factorial(UU + UD))/10^20)* ((factorial(DD + DU))/10^20)* ((factorial(UU + DU))/10^20)* ((factorial(DD + UD))/10^20))/(((factorial(UU))/10^16)*((factorial(DD))/10^16)*((factorial(UD))/10^16)*((factorial(DU))/10^16)*((factorial(UU+DD+UD+DU))/10^16));

%% table for IEMs
Biomarkers = unique(R22.Swagatika_IEM_biomarker(:,2));
IEMsSS = unique(R22.Swagatika_IEM_biomarker(:,1));
clear IEMMatrixR22 IEMxR22
cntR=2;
cntC=2;
cnt = 1;
IEMMatrixR22(1,1)={''};
for i =1  : size(R22.Swagatika_IEM_biomarker,1)
    Change = str2num(R22.Swagatika_IEM_biomarker{i,3});
    if Change ~= 0
        IEMxR22(cnt,1)=R22.Swagatika_IEM_biomarker(i,1);
        cnt = cnt+1;
    end
end

IEMxR22 = unique(IEMxR22);

% this matrix contains  predicted and reported biomarkers
IEMMatrixR22(2:length(IEMxR22)+1,1)=IEMxR22;
for i = 1 : size(R22.Swagatika_IEM_biomarker,1)
    Change = str2num(R22.Swagatika_IEM_biomarker{i,3});
    R22Change = R22.Swagatika_IEM_biomarker{i,5};
    % if R2Change ~= 0  || ~isempty(strmatch(R2.Swagatika_IEM_biomarker(i),IEMMatrix(1,:),'exact')) || ~isempty(strmatch(R2.Swagatika_IEM_biomarker(i),IEMMatrix(:,1),'exact'))
    C=strmatch(R22.Swagatika_IEM_biomarker(i,2),IEMMatrixR22(1,:),'exact');
    R=strmatch(R22.Swagatika_IEM_biomarker(i,1),IEMMatrixR22(:,1),'exact');
    if ~isempty(R)
        if isempty(C)
            IEMMatrixR22(1,cntC)= R22.Swagatika_IEM_biomarker(i,2);
            C= cntC;cntC = cntC+1;
        end
        if R22Change ==2
            R22Change =1;
        elseif R22Change== -2
            R22Change=-1;
        end
        if Change ==R22Change
            if Change==-1
                IEMMatrixR22(R,C)={'DD'};
            elseif  Change==1
                IEMMatrixR22(R,C)={'UU'};
            end
        elseif  Change==-1 && R22Change == 1
            IEMMatrixR22(R,C)={'DU'};
        elseif  Change==1 && R22Change == -1
            IEMMatrixR22(R,C)={'UD'};
        elseif  Change==1 && R22Change == 0
            IEMMatrixR22(R,C)={'U0'};
        elseif  Change==-1 && R22Change == 0
            IEMMatrixR22(R,C)={'D0'};
        end
    end
end

G = R22.IEMGenes(ismember(R22.IEMNames,IEMxR22));
G = strrep(strcat('G_',G), ':','_');
I =  find(ismember(R22.IEMsList(:,1),G));
x=1;
for i = 1 : length(I)
    C=strmatch(R22.IEMsList(I(i),2),IEMMatrixR22(1,:),'exact');
    if ~isempty(C)
        TMP(x,1:3)=R22.IEMsList(I(i),1:3);

        R22Change = R22.IEMsList{I(i),3};
        X=regexprep(R22.IEMsList(I(i),1),'G_','');
        Y=strmatch(strrep(X,'_',':'),R22.IEMGenes,'exact');
        TMP(x,4)=R22.IEMNames(Y(1));
        x=x+1;
        R=strmatch(R22.IEMNames(Y(1)),IEMMatrixR22(:,1),'exact');
        if ~isempty(R) && isempty(IEMMatrixR22{R,C})
            if R22Change <0
                IEMMatrixR22(R,C)={'0D'};
            elseif R22Change >0
                IEMMatrixR22(R,C)={'0U'};
            end
        end
    end
end

% this matrix contains only those biomarkers that have been reported in the
% literature
clear IEMMatrix2R22
IEMMatrix2R22(1,1)={''};
IEMMatrix2R22(2:length(IEMxR22)+1,1)=IEMxR22;
cntR=2;
cntC=2;
for i = 1 : size(R22.Swagatika_IEM_biomarker,1)
    Change = str2num(R22.Swagatika_IEM_biomarker{i,3});
    R22Change = R22.Swagatika_IEM_biomarker{i,5};
    % if R2Change ~= 0  || ~isempty(strmatch(R2.Swagatika_IEM_biomarker(i),IEMMatrix2(1,:),'exact')) || ~isempty(strmatch(R2.Swagatika_IEM_biomarker(i),IEMMatrix2(:,1),'exact'))
    C=strmatch(R22.Swagatika_IEM_biomarker(i,2),IEMMatrix2R22(1,:),'exact');
    R=strmatch(R22.Swagatika_IEM_biomarker(i,1),IEMMatrix2R22(:,1),'exact');
    if ~isempty(R)
        if  Change==1 && R22Change == 0 ||Change==-1 && R22Change == 0
        else
            if isempty(C)
                IEMMatrix2R22(1,cntC)= R22.Swagatika_IEM_biomarker(i,2);
                C= cntC;cntC = cntC+1;
            end
            if R22Change ==2
                R22Change =1;
            elseif R22Change== -2
                R22Change=-1;
            end
            if Change ==R22Change
                if Change==-1
                    IEMMatrix2R22(R,C)={'DD'};
                elseif  Change==1
                    IEMMatrix2R22(R,C)={'UU'};
                end
            elseif  Change==-1 && R22Change == 1
                IEMMatrix2R22(R,C)={'DU'};
            elseif  Change==1 && R22Change == -1
                IEMMatrix2R22(R,C)={'UD'};
            end
        end
    end
end

save('ResultsIEManalysisR22_flavo.mat');
% clear unneccessary variables.
clear G I Genes IEMs* R NBB List IEMs R2BD* U* i b* a* Y X c cnt* cons* dele* has* tok rem modelR2* solv* upreg* x tmp Ex* FR* D* Bioma* C* j TMP R2B* R2G* R2C* IEMx