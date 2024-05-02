

ST_Data_Integration <- function(visium_dir = visium_dir, method = "FeatureAnchoring", conda_env = "/miniconda3/envs/scvi",epochs = 100) {
    if (!require("batchelor", character.only = TRUE)) {
        BiocManager::install("batchelor")
    }

    if (!require("Seurat", character.only = TRUE)) {
        install.packages("Seurat")
    }

    if (!require("future", character.only = TRUE)) {
        install.packages("future")
    }
    if (!require("SingleCellExperiment", character.only = TRUE)) {
     BiocManager::install("SingleCellExperiment")
    }
    if (!require("reticulate", character.only = TRUE)) {
        install.packages("reticulate")
    }

    file_names <- list.files(path = visium_dir)
    section_ids <- c(file_names)

    for (i in section_ids) {
        path <- paste0(visium_dir, i, "/outs")
        im <- Read10X_Image(
            paste0(path, "/spatial"),
            image.name = "tissue_lowres_image.png",
            filter.matrix = TRUE,
        )
        eval(parse(text = paste0(i, " <- Load10X_Spatial(data.dir = path, image = im)")))
    }
    # create a list of the original data that we loaded to start with
    st.list <- mget(section_ids)
    for (idx in seq_along(st.list)) {
        st.list[[idx]] <- st.list[[idx]][, unname(which(colSums(GetAssayData(st.list[[idx]])) != 0))] # remove all zero spots
        st.list[[idx]] <- subset(st.list[[idx]], nFeature_Spatial > 100 & nCount_Spatial > 0) # filter out low quality spots
        st.list[[idx]][["orig.ident"]] <- section_ids[idx] # add orig.ident to label each slice
        st.list[[idx]] <- UpdateSeuratObject(st.list[[idx]])
    }

    if (!method %in% c("FeatureAnchoring", "SCVI", "FastMNN_variable","FastMNN_all", "CCA", "Harmony")) {
        print(" please select method  from the following options: 'FeatureAnchoring','SCVI','FastMNN_variable','FastMNN_all','CCA','Harmony'")
    }
    ## Integrate multiple slice
    # run SCT on both datasets
    st.list <- lapply(st.list, SCTransform, vst.flavor = "v2", assay = "Spatial", method = "poisson")

    if (method == "FastNMF_all" | method == "FastNMF_variable" ) {
        ## FastNMF integration
        # BiocManager::install("batchelor")
        merged_spe <- Reduce(merge, st.list)
        merged_spe <- SetIdent(merged_spe, value = "orig.ident")
        layers <- Layers(merged_spe, search = "data")
        if (method == "FastNMF_all") {
            cm_features <- Reduce(intersect, lapply(st.list, Features))
        } else if (method == "FastNMF_variable" ) {
            cm_features <- Reduce(intersect, lapply(st.list, VariableFeatures))
        }
        counts <- as.matrix(LayerData(merged_spe, assay = "SCT", layer = "data")[cm_features, colnames(LayerData(merged_spe, layer = layers))])
        combined <- SingleCellExperiment(list(counts = counts), metadata = data.frame(merged_spe[[]]))
        counts <- assay(combined, "counts")
        libsizes <- colSums(counts)
        size.factors <- libsizes / mean(libsizes)
        logcounts(combined) <- log2(t(t(counts) / size.factors) + 1)
        metadata <- data.frame(merged_spe[[]])
        out <- batchelor::fastMNN(combined,
            batch = metadata$orig.ident,
            auto.merge = TRUE,
            assay.type = "counts"
        )
        merged_spe[["integrated_NMF"]] <- CreateDimReducObject(
            embeddings = SingleCellExperiment::reducedDim(x = out),
            loadings = as.matrix(SingleCellExperiment::rowData(x = out)),
            assay = DefaultAssay(object = merged_spe),
            key = "NMF_",
        )
        integ_spe <- FindNeighbors(merged_spe, dims = 1:50, reduction = "integrated_NMF")
        integ_spe <- FindClusters(integ_spe, resolution = 2, cluster.name = "integrated_clusters")
        integ_spe <- RunUMAP(integ_spe, dims = 1:50, reduction = "integrated_NMF", reduction.name = "umap.NMF")
        # DimPlot(integ_spe, reduction = "umap.NMF", group.by = c("integrated_clusters", "orig.ident"), pt.size = 3)
    } else if (method == "Harmony" | method == "CCA") {
        # integrated AS layers
        # Example to merge more than two Seurat objects
        merged_spe <- Reduce(merge, st.list)
        merged_spe <- SetIdent(merged_spe, value = "orig.ident")

        DefaultAssay(merged_spe) <- "SCT"
        merged_spe <- NormalizeData(merged_spe)
        VariableFeatures(merged_spe) <- unique(unlist(lapply(st.list, VariableFeatures)))
        merged_spe <- ScaleData(merged_spe)
        merged_spe <- RunPCA(merged_spe)
        if (method == "CCA") {
            integ_spe <- IntegrateLayers(
                object = merged_spe, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                verbose = FALSE
            )
        }
        if (method == "Harmony") {
            integ_spe <- IntegrateLayers(
                object = merged_spe,
                method = HarmonyIntegration, assay = "SCT",
                orig.reduction = "pca", new.reduction = "harmony",
                verbose = TRUE
            )
        }
        integ_spe <- FindNeighbors(integ_spe, dims = 1:30, reduction = "harmony")
        integ_spe <- FindClusters(integ_spe, resolution = 2, cluster.name = "integrated_clusters")
        integ_spe <- RunUMAP(integ_spe, dims = 1:30, reduction = "harmony", reduction.name = "umap.integrated")
    } else if (method == "FeatureAnchoring") {
        # need to set maxSize for PrepSCTIntegration to work
        options(future.globals.maxSize = 2000 * 1024^2) # set allowed size to 2K MiB
        cm_var_features <- intersect(Reduce(intersect, lapply(st.list, Features)), Reduce(intersect, lapply(st.list, VariableFeatures)))
        st.features <- SelectIntegrationFeatures(object.list = st.list, nfeatures = length(cm_var_features))
        st.list <- PrepSCTIntegration(object.list = st.list, anchor.features = length(cm_var_features), verbose = FALSE)
        int.anchors <- FindIntegrationAnchors(object.list = st.list, normalization.method = "SCT", verbose = FALSE, anchor.features = st.features)
        integ_spe <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)
        integ_spe <- RunPCA(integ_spe, verbose = FALSE)
        integ_spe <- FindNeighbors(integ_spe, dims = 1:30)
        integ_spe <- FindClusters(integ_spe, verbose = FALSE)
        integ_spe <- RunUMAP(integ_spe, dims = 1:30)
        # DimPlot(integ_spe, reduction = "umap", group.by = c("ident", "orig.ident"))
        # SpatialDimPlot(integ_spe, pt.size.factor = 3)
        # Loadings(integ_spe, reduction = "pca")
        # Embeddings(integ_spe, reduction = "pca")
    } else if (method == "SCVI") {
        ## SCVI integration
        if (conda_env == "") {
            print("Please specify conda env")
        }
        # import python methods from specified conda env
        reticulate::use_condaenv(conda_env, required = TRUE)
        sc <- reticulate::import("scanpy", convert = FALSE)
        scvi <- reticulate::import("scvi", convert = FALSE)
        anndata <- reticulate::import("anndata", convert = FALSE)
        scipy <- reticulate::import("scipy", convert = FALSE)
        merged_spe <- Reduce(merge, st.list)
        # build a meta.data-style data.frame indicating the batch for each cell
        batches <- merged_spe[[]] |>
            dplyr::select(orig.ident) |>
            dplyr::rename(batch = orig.ident)
        # setup an `AnnData` python instance
        adata <- sc$AnnData(
            X = scipy$sparse$csr_matrix(Matrix::t(LayerData(merged_spe, layer = "counts"))),
            obs = batches,
            var = data.frame(row.names = Features(merged_spe, layer = "counts")),
        )
        scvi$model$SCVI$setup_anndata(adata, batch_key = "batch")

        # initialize and train the model
        model <- scvi$model$SCVI(
            adata = adata, n_latent = as.integer(10), n_layers = as.integer(1), gene_likelihood = "nb"
        )
        model$train(max_epochs = as.integer(epochs))
        # extract the latent representation of the merged data
        latent <- model$get_latent_representation()
        latent <- as.matrix(latent)
        # pull the cell identifiers back out of the `AnnData` instance
        rownames(latent) <- reticulate::py_to_r(adata$obs$index$values)
        colnames(latent) <- paste0("integrated.scvi", "_", 1:ncol(latent))

        # build a `DimReduc` instance
        merged_spe[["integrated_scvi"]] <- CreateDimReducObject(embeddings = latent, key = "integrated_scvi", assay = DefaultAssay(merged_spe))
        integ_spe <- FindNeighbors(merged_spe, dims = 1:ncol(latent), reduction = "integrated_scvi")
        integ_spe <- FindClusters(integ_spe, resolution = 2, cluster.name = "integrated_clusters")
        integ_spe <- RunUMAP(integ_spe, dims = 1:ncol(latent), reduction = "integrated_scvi", reduction.name = "umap.integrated")
        # DimPlot(integ_spe, reduction = "umap.integrated", group.by = c("integrated_clusters", "orig.ident"), pt.size = 3)
    }

}




