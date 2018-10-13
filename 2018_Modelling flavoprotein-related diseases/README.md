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

Models have been sampled using optGpSampler.

<b>Fig. 4. and Fig. S4.</b> Coupling of FAD-related reactions to FAD-biosynthesis enabled the new model to respond to low cofactor availability.

To introduce cofactor dependency in the FAD model use: cofactorLimitation.m . By changing the cofactor coefficience in the line 43 we can test sensitivity of the metabolism to the cofactor coefficience value. To calculate the average flux through the flavoprotein-dependent reaction run the cofactorLimitationTest.m (Fig. 4B, Fig. S4). To calculate the sampled average flux through the flavoprotein-dependent reaction run the cofactorSampling.m (Fig. 4C). optGpSampler is required to run the cofactorSampling.m script. 

<b>Fig. 5.</b> Changes in the ATP yield from different carbon sources in flavoprotein-related diseases predict metabolic adaptations in energy metabolism.

<b>Table. 2, Table. S3, and Fig. S3</b> Metabolic biomarkers for flavoprotein-related diseases.

