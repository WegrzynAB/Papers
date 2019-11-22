%// Your code
samples_C = loadSamples('control_fibro', 10, 1000, 0);
samples_R = loadSamples('refsum_fibro',  10, 1000, 0);
samples_C_phyt = loadSamples('control_fibro_phyt', 10, 1000, 0);
samples_R_phyt = loadSamples('refsum_fibro_phyt', 10, 1000, 0);

S_fibroblast_C = samples_C(exIdx,:);
S_fibroblast_R = samples_R(exIdx,:);
S_fibroblast_C_phyt = samples_C_phyt(exIdx,:);
S_fibroblast_R_phyt = samples_R_phyt(exIdx,:);

all_samples = [S_fibroblast_C'; S_fibroblast_R'; S_fibroblast_C_phyt'; S_fibroblast_R_phyt'];

[coeff,score,latent,tsquared,explained] = pca(all_samples);    
score = score';

%// Group 1: CTRL
group1 = score(1:2,1:10000);
%// Group 2: RD
group2 = score(1:2,10001:20000);
%// Group 3: CTRL + phyt
group3 = score(1:2,20001:30000);
%// Group 4: RD + phyt
group4 = score(1:2,30001:40000); 

%// Plot as separate colours
hold on
scatter(group1(1,:), group1(2,:), 'MarkerEdgeColor', [0.7 0.7 0.7], 'MarkerFaceColor', [0.7 0.7 0.7], 'Marker', '.', 'MarkerFaceAlpha', .5);
scatter(group2(1,:), group2(2,:), 'MarkerEdgeColor',[0 0 0],'MarkerFaceColor', [0 0 0], 'Marker', '.', 'MarkerFaceAlpha', .5); 
scatter(group3(1,:), group3(2,:), 'MarkerEdgeColor', [0.9 0.9 0.9], 'MarkerFaceColor',[0.9 0.9 0.9],'Marker','.', 'MarkerFaceAlpha', .5); 
scatter(group4(1,:), group4(2,:), 'MarkerEdgeColor',[0.4 0.4 0.4], 'MarkerFaceColor',[0.4 0.4 0.4],'Marker', '.', 'MarkerFaceAlpha', .5);   
hold off
