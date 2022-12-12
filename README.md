# Trajectory_Analysis
1) makeSeurat.R creates the Seurat object and performs dimensionality reduction.
2) makeSlingshot.R finds trajectories based on dimensionality reduction.
3) getMetadata.R creates U matrix to correct for covariates in fitGAM.R.
4) makeJob.R generates joblist.txt for dead simple queue based on number of genes 
   and desired batch size.
5) fitGAM.R is to be run by the dead simple queue (you should not have to open this file).
6) combineGAM.R will combine the output from dead simple queue process.
7) dsAnalysis.R will get DE genes across lineages.
8) compareLins.R will compare computed lineages from each dataset (will need to edit based on additional
   datasets to be used). 
