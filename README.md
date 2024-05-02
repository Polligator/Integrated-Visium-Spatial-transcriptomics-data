---
title: "Readme"
output: github_document
---


A simple function used to integrate 10x Genomics Visium Spatial transcriptomics data, the integration is based on Seuart frame work.
integration methods includes: scVI,FastMNN,CCA,Harmony and FeatureAnchoring.

To use this function, you need to have the following R packages installed: Seurat, batchelor, reticulate, and the following python packages: scanpy, anndata, scvi-tools.
All you need to do is to provide the path to the directory that contains all the spaceranger output data, the function will automatically read the data and integrate them.
here is the structure of the directory should look like for this function to work:

```
Result
├── Sample1
│   ├── outs
│   │   ├── analysis
│   │   ├── filtered_feature_bc_matrix
│   │   ├── raw_feature_bc_matrix
│   │   ├── spatial
│   │       ├── tissue_positions.csv
│   │       ├── scalefactors_json.json
│   │       ├── image.tif
│   │       ├── tissue_hires_image.png
│   │       ├── tissue_lowres_image.png
│   │       ├── aligned_fiducials.jpg
│   │       ├── detected_tissue_image.jpg
│   │       ├── spatial_enrichment.csv
├── Sample2
├── Sample3
│--- Sample4
│--- Sample5
│--- Sample6
│--- Sample7
│--- Sample8
│
```
oragnize your data in this way, and provide the path to the Result directory, the function will automatically read the data and integrate them.


## scVI: 
To use this method, you need to have a working conda environment set up for the scVI, please set conda_env=Path/to/your/env, for example conda_env="/home/USER/miniconda3/envs/scvi"

To load this function and run in your R terminal:

```{r setup, include=FALSE}
source('https://raw.githubusercontent.com/Polligator/Integrated-10x-Genomics-Visium-Spatial-transcriptomics-data/main/integration.r')
inetgrated_ST<-ST_Data_Integration(visium_dir = visium_dir, method = "SCVI", conda_env = "/miniconda3/envs/scvi",epochs = 100) 
```
visium_dir is the path to your directory, which should contain all the individual spaceranger data folder.

To Visualize:
```{r setup, include=FALSE}
DimPlot(inetgrated_ST, reduction = "umap.integrated", group.by = c("integrated_clusters", "orig.ident"), pt.size = 3)
```

## FastMNN: 
This method integrate ST data using batchelor packages : https://bioconductor.org/packages/release/bioc/html/batchelor.html, there are two options for this method: FastMNN_all and FastNMF_variable. FastMNN_all using all common genes that present in all you samples to integrate the data, FastNMF_variable only use the vaibale features that identified by Seurat function VariableFeatures() to integrate

```{r setup, include=FALSE}
inetgrated_ST <- ST_Data_Integration(visium_dir = visium_dir, method = "FastNMF_variable")
inetgrated_ST <- ST_Data_Integration(visium_dir = visium_dir, method = "FastNMF_all")

```
To Visualize:
```{r setup, include=FALSE}
DimPlot(inetgrated_ST, reduction = "umap.NMF", group.by = c("integrated_clusters", "orig.ident"), pt.size = 3)

```


## CCA or Harmony

```{r setup, include=FALSE}
inetgrated_ST <- ST_Data_Integration(visium_dir = visium_dir, method = "CCA")
inetgrated_ST <- ST_Data_Integration(visium_dir = visium_dir, method = "Harmony")
``` 
To Visualize:
```{r setup, include=FALSE}
DimPlot(inetgrated_ST, reduction = "umap.integrated", group.by = c("integrated_clusters", "orig.ident"), pt.size = 3)
```

## FeatureAnchoring
This method is based on the Seurat package, it uses the feature anchoring function to integrate the data, the function will automatically select the maximiun number of variable features to integrate the data.
Esstentially, it is a wrapper function for the Seurat function FindIntegrationAnchors and IntegrateData.
```{r setup, include=FALSE}
int.anchors <- FindIntegrationAnchors(object.list = st.list, normalization.method = "SCT", verbose = FALSE, anchor.features = st.features)
brain.integrated <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)
```
To use this method
```{r setup, include=FALSE}
inetgrated_ST <- ST_Data_Integration(visium_dir = visium_dir, method = "FeatureAnchoring")
```

To Visualize:
```{r setup, include=FALSE}
DimPlot(inetgrated_ST, reduction = "umap", group.by = c("ident", "orig.ident"))

```
