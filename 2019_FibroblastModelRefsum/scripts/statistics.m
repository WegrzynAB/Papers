models = {'fibroblast_C_phyt', 'fibroblast_R_phyt'}; %, 'hexoK'
exIdx = findExcRxns(fibroblast_C,'true');
exRxns = fibroblast_C_phyt.rxns(exIdx);
S_fibroblast_C_phyt = samples_C_phyt(exIdx,:);
S_fibroblast_R_phyt = samples_R_phyt(exIdx,:);
S_fibroblast_C = samples_C(exIdx,:);
S_fibroblast_R = samples_R(exIdx,:);
%, 'hexoK'
%statistical test for differences between samples
mean_C = mean(S_fibroblast_C_phyt,2);
mean_R = mean(S_fibroblast_R_phyt,2);
mean_FC = mean_R./mean_C;

median_C = median(S_fibroblast_C_phyt,2);
median_Rp = median(S_fibroblast_R_phyt,2);
mean_Rc = mean(S_fibroblast_R,2);

for i = 1:length(exRxns)
    if mean_FC(i) >= 0
        changeDirect(i) = 0;
        ksTest_C(i) = kstest(S_fibroblast_C_phyt(i,:));
        ksTest_R(i) = kstest(S_fibroblast_R_phyt(i,:));
        [~, ksTest2_p(i)] = kstest2(S_fibroblast_C_phyt(i,:),S_fibroblast_R_phyt(i,:));
        if ksTest_C(i) == 0 && ksTest_R == 0
            changeDist(i) = 0;
            Ftest_h(i) = vartest2(S_fibroblast_C_phyt(i,:),S_fibroblast_R_phyt(i,:));
            if Ftest_h(i) == 0
                [ttest_h(i), ttest_p(i)] = ttest2(S_fibroblast_C_phyt(i,:),S_fibroblast_R_phyt(i,:));
            else
                [ttest_h(i), ttest_p(i)] = ttest2(S_fibroblast_C_phyt(i,:),S_fibroblast_R_phyt(i,:), 'Vartype','unequal');
            end
        elseif ksTest_C(i) == ksTest_R(i) %non-normal distribution in both
            changeDist(i) = 0;
            [wilcox_p(i),wilcox_h(i)] = ranksum(S_fibroblast_C_phyt(i,:),S_fibroblast_R_phyt(i,:));
        else %change between normal and non-normal distribution
            changeDist(i) = 1;
            [wilcox_p(i),wilcox_h(i)] = ranksum(S_fibroblast_C_phyt(i,:),S_fibroblast_R_phyt(i,:));
        end
    else
        changeDirect(i) = 1;
    end
end

[FDR, h] = bonf_holm(wilcox_p);
[FDR2, h2] = bonf_holm(ksTest2_p);

wilcox_p = wilcox_p';
FDR_wilcox = FDR';
ksTest2_p = ksTest2_p';
FDR_ks = FDR2';
log_FC = log2(mean_FC);


table_res = table(exRxns, mean_C, mean_R, mean_FC, log_FC, wilcox_p, FDR_wilcox, ksTest2_p, FDR_ks,'RowNames', exRxns);
writetable(table_res, 'stats_table_3.txt')

Phyt = {'EX_phyt[e]';'EX_prist[e]'};

CTRL = [];
RD = [];
metNames = {};
for i=1:length(Phyt)
    ID = findRxnIDs(fibroblast_C_phyt, Phyt(i));
    CTRL = [CTRL samples_C_phyt(ID,:)'];
    RD = [RD samples_R_phyt(ID,:)'];
    metNames(i) = fibroblast_C_phyt.metNames(findMetIDs(fibroblast_C_phyt, findMetsFromRxns(fibroblast_C_phyt, Phyt(i))));
end

metNames = metNames';

distributionPlot(CTRL,'histOri','left','color',[0.85 0.85 0.85],'widthDiv',[2 1],'showMM',0,'histOpt',2, 'xNames', metNames)
distributionPlot(RD,'histOri','right','color',[0 0 0],'widthDiv',[2 2],'showMM',0, 'histOpt',2)

MAA = {'EX_3MAA[e]';'EX_dmhptcrn[e]'};

CTRL = [];
RD = [];
metNames = {};
for i=1:length(MAA)
    ID = findRxnIDs(fibroblast_C_phyt, MAA(i));
    CTRL = [CTRL samples_C_phyt(ID,:)'];
    RD = [RD samples_R_phyt(ID,:)'];
    metNames(i) = fibroblast_C_phyt.metNames(findMetIDs(fibroblast_C_phyt, findMetsFromRxns(fibroblast_C_phyt, MAA(i))));
end

metNames = metNames';

distributionPlot(CTRL,'histOri','left','color',[0.85 0.85 0.85],'widthDiv',[2 1],'showMM',0,'histOpt',2, 'xNames', metNames)
distributionPlot(RD,'histOri','right','color',[0 0 0],'widthDiv',[2 2],'showMM',0, 'histOpt',2)

table_sig = table_res.exRxns(table_res.FDR_wilcox < 0.05 & abs(table_res.log_FC) > 1.3);

for i=1:length(table_res.exRxns)
    names(i,1) = fibroblast_C.metNames(findMetIDs(fibroblast_C,findMetsFromRxns(fibroblast_C, table_res.exRxns(i))));
    names(i,2) = table_res.exRxns(i);
    names(i,3) = findMetsFromRxns(fibroblast_C, table_res.exRxns(i));
end


CTRL = [];
RD = [];
for i=1:length(table_sig)
    ID = findRxnIDs(fibroblast_C_phyt, table_sig(i));
    CTRL = [CTRL samples_C_phyt(ID,:)'];
    RD = [RD samples_R_phyt(ID,:)'];
    metNames(i) = fibroblast_C_phyt.metNames(findMetIDs(fibroblast_C_phyt, findMetsFromRxns(fibroblast_C_phyt, table_sig(i))));
end

metNames = metNames';

distributionPlot(CTRL,'histOri','left','color',[0.85 0.85 0.85],'widthDiv',[2 1],'showMM',0,'histOpt',2, 'xNames', metNames)
distributionPlot(RD,'histOri','right','color',[0 0 0],'widthDiv',[2 2],'showMM',0, 'histOpt',2)

tablesig_2 = {'EX_fe3[e]';'EX_lac_D[e]';'EX_ala_L[e]';'EX_mercplaccys[e]'};
CTRL = [];
RD = [];
metNames = {};
for i=1:length(tablesig_2)
    ID = findRxnIDs(fibroblast_C_phyt, tablesig_2(i));
    CTRL = [CTRL samples_C_phyt(ID,:)'];
    RD = [RD samples_R_phyt(ID,:)'];
    metNames(i) = fibroblast_C_phyt.metNames(findMetIDs(fibroblast_C_phyt, findMetsFromRxns(fibroblast_C_phyt, tablesig_2(i))));
end

metNames = metNames';

distributionPlot(CTRL,'histOri','left','color',[0.85 0.85 0.85],'widthDiv',[2 1],'showMM',0,'histOpt',2, 'xNames', metNames)
distributionPlot(RD,'histOri','right','color',[0 0 0],'widthDiv',[2 2],'showMM',0, 'histOpt',2)

C6 = {'EX_caproic[e]';'EX_citr_L[e]'}; %;'EX_3thexddcoacrn[e]'

CTRL = [];
RD = [];
metNames = {};
for i=1:length(C6)
    ID = findRxnIDs(fibroblast_C_phyt, C6(i));
    CTRL = [CTRL samples_C_phyt(ID,:)'];
    RD = [RD samples_R_phyt(ID,:)'];
    metNames(i) = fibroblast_C_phyt.metNames(findMetIDs(fibroblast_C_phyt, findMetsFromRxns(fibroblast_C_phyt, C6(i))));
end

metNames = metNames';

distributionPlot(CTRL,'histOri','left','color',[0.85 0.85 0.85],'widthDiv',[2 1],'showMM',0,'histOpt',2, 'xNames', metNames)
distributionPlot(RD,'histOri','right','color',[0 0 0],'widthDiv',[2 2],'showMM',0, 'histOpt',2)

tablesig_2 = {'EX_CE1557[e]';'EX_argcysser[e]';'EX_gluilelys[e]';'EX_proproarg[e]';'EX_trpthrglu[e]';'EX_tyrthr[e]';'EX_valhisasn[e]'};
CTRL = [];
RD = [];
metNames = {};
for i=1:length(tablesig_2)
    ID = findRxnIDs(fibroblast_C_phyt, tablesig_2(i));
    CTRL = [CTRL samples_C_phyt(ID,:)'];
    RD = [RD samples_R_phyt(ID,:)'];
    metNames(i) = fibroblast_C_phyt.metNames(findMetIDs(fibroblast_C_phyt, findMetsFromRxns(fibroblast_C_phyt, tablesig_2(i))));
end

metNames = metNames';

distributionPlot(CTRL,'histOri','left','color',[0.85 0.85 0.85],'widthDiv',[2 1],'showMM',0,'histOpt',2, 'xNames', metNames)
distributionPlot(RD,'histOri','right','color',[0 0 0],'widthDiv',[2 2],'showMM',0, 'histOpt',2)
xtickangle(45)

