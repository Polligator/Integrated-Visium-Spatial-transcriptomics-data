A simple function used to integrate 10x Genomics Visium Spatial transcriptomics data, the integration is based on Seuart frame work, integration methods includes: scVI,FastMNN,CCA,Harmony and Feature_Anchoring.

scVI: to use this method, you need to have a working conda environment set up for the scVI, please set conda_env=Path/to/your/env, for example conda_env="/miniconda3/envs/scvi"

To load this function and run in your R terminal:

source('https://raw.githubusercontent.com/Polligator/Integrated-10x-Genomics-Visium-Spatial-transcriptomics-data/main/integration.r')
inetgrated_ST<-ST_Data_Integration(visium_dir = visium_dir, method = "SCVI", conda_env = "/miniconda3/envs/scvi")
visium_dir is the path to your directory, which should contain all the individual spaceranger data folder.

FastMNN: this method integrate ST data using batchelor packages : https://bioconductor.org/packages/release/bioc/html/batchelor.html, there are two options for this method: FastMNN_all and FastMNN_vaiable. FastMNN_all using all common genes that present in all you samples to integrate the data, FastMNN_vaiable only use the vaibale features that identified by Seurat function VariableFeatures() to integrate




