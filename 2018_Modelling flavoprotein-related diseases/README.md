The scripts in this folder are relative to:

<b>Cofactors revisited – predicting the impact of flavoprotein-related diseases on a genome scale.</b>

Agnieszka B. Wegrzyn, Sarah Stolle, Rienk A. Rienksma, Vítor A. P. Martins dos Santos, Barbara M. Bakker$, and Maria Suarez-Diez$

$These authors contributed equally. 

Available at: https://www.biorxiv.org/content/early/2018/08/15/343624

<b>Additional packages required</b>

To reproduce the results from our study you need COBRA Toolbox https://opencobra.github.io/cobratoolbox/stable/, and optGpSampler http://cs.ru.nl/~wmegchel/optGpSampler/

<b>Fig. 2. and Fig. S1.</b> Mapping of the flavoproteome on the models before (Recon 2.2) and after curation (Recon 2.2_FAD).

To reproduce the data please run the Fig2_flavoproteomeMapping.m script. Table S1 contains information about the flavoprotein genes and names required to run the script (loaded withing the script). 

<b>Fig. 3., Fig. S2., and Fig S5.</b> New models can correctly simulate the physiology of MADD.

For Recon2.2 MADD simulations use Fig3_MADDsamplingRecon22.m script (Fig.3 and Fig. S2). For Recon3D MADD simulations use FigS5_MADDsamplingRecon3D script (Fig. S5). For both scripts 'glpk' or 'gurobi' solver is required together with the optGpSampler.

<b>Fig. 4. and Fig. S4.</b> Coupling of FAD-related reactions to FAD-biosynthesis enabled the new model to respond to low cofactor availability.

To introduce cofactor dependency in the FAD model use: Fig4_cofactorLimitation.m . By changing the cofactor coefficience in the line 43 we can test sensitivity of the metabolism to the cofactor coefficience value. To calculate the average flux through the flavoprotein-dependent reaction run the Fig4_cofactorLimitationTest.m (Fig. 4B, Fig. S4). To calculate the sampled average flux through the flavoprotein-dependent reaction run the Fig4c_cofactorSampling.m (Fig. 4C). optGpSampler is required to run the Fig4c_cofactorSampling.m script. 

<b>Fig. 5.</b> Changes in the ATP yield from different carbon sources in flavoprotein-related diseases predict metabolic adaptations in energy metabolism.

<b>Table. 2, Table. S3, and Fig. S3</b> Metabolic biomarkers for flavoprotein-related diseases.

For Recon2.04 and Recon3D we used the original method and files linked to the Recon2 paper (https://github.com/opencobra/COBRA.papers , 2013_Recon2 folder). For Recon2.2 we have slightly updated the scripts to work with HGNC gene identifiers, see T2_IEMAnalysisRecon22.m script. Modified script requires two files: IEM_biomarkerListR22.mat  and IEM_compendiumR22.mat that contain information about genes and biomarkers asociated with the diseases. 
