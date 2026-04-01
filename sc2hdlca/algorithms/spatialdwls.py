def spatialdwls(
    adata_spatial, 
    adata_ref, 
    labels_key,
    cluster_key = "leiden",
    n_cell=50,
    r_lib_path=None, 
    python_path=None,
    results_path="./spatialdwls_results", 
): 
    
    '''
    Run SpatialDWLS

    Parameters 
    -----------

    adata_spatial : AnnData
        AnnData of the spatial data , filtered by highly variable features. Feature space needs to be the same as the one of adata_ref.
        Normalized counts are expected to be saved in .X. Names of observations and features are expected as row names of adata_spatial.obs and adata_spatial.var. 
    adata_ref : AnnData 
        AnnData of the reference data, filtered by highly variable features. Feature space needs to be the same as the one of adata_spatial.
         Normalized counts are expected to be saved in .X. Names of observations and features are expected as row names of adata_ref.obs and adata_ref.var. 
    labels_key : str
        Cell type key in adata_ref.obs for label information.
    cluster_key : str
        Cluster key in adata_spatial.obs for cluster information.
    n_cell : int
        Number of cells per spot.
    r_lib_path : str
        Path to R library.  
    python_path : str
        Path to R library.  
    results_path : str
        Path to save estimated cell type abundances to. 
    
        
    Returns
    --------

    - Saves estimated proportions as csv-file to results_path.
    '''

    from rpy2 import robjects
    import anndata2ri
    import scanpy as sc
    
    adata_spatial_copy = adata_spatial.copy()
    adata_ref_copy = adata_ref.copy()

    if r_lib_path is not None:
        robjects.r.assign("r_path", r_lib_path)
        robjects.r(".libPaths(r_path)")

    # Load reference into R and create signature matrix: 
    scexp_sc = anndata2ri._py2r.py2rpy_anndata(adata_ref_copy)
    robjects.r.assign("scexp_sc", scexp_sc)
    robjects.r.assign("labels_key", labels_key)
    robjects.r("""
        library(Giotto)
        library(data.table)
        library(SingleCellExperiment)
    """)

    # =========================
    # Reference (scRNA-seq)
    # =========================
    print("Loading reference into R and creating signature matrix")

    scexp_sc = anndata2ri._py2r.py2rpy_anndata(adata_ref_copy)
    robjects.r.assign("scexp_sc", scexp_sc)
    robjects.r.assign("labels_key", labels_key)
    robjects.r.assign("python_path", python_path)
    robjects.r("""
                # instructions
                my_python_path = python_path
                instrs = createGiottoInstructions(python_path = my_python_path)
        
                # expression
                counts_sc = assay(scexp_sc, "X")
        
                # metadata
                meta_sc = as.data.table(colData(scexp_sc))
                meta_sc$cell_ID = rownames(colData(scexp_sc))
                meta_sc$new_celltype = meta_sc[[labels_key]]
        
                var_sc = as.data.table(rowData(scexp_sc))
                var_sc$feat_ID = rownames(rowData(scexp_sc))
        
                # Giotto object
                sc_obj = createGiottoObject(
                    raw_exprs     = counts_sc,
                    instructions  = instrs,
                    cell_metadata = meta_sc,
                    gene_metadata = var_sc
                )
        
                # preprocess
                sc_obj = filterGiotto(
                    gobject = sc_obj,
                    expression_threshold = 0.1,
                    gene_det_in_min_cells = 10,
                    min_det_genes_per_cell = 10,
                    expression_values = c("raw"),
                    verbose = TRUE
                )
        
                sc_obj = normalizeGiotto(
                    gobject = sc_obj,
                    scalefactor = 6000,
                    verbose = TRUE
                )
        
                sc_obj = addStatistics(sc_obj)
        
                # gini markers
                gini_markers = findMarkers_one_vs_all(
                    gobject = sc_obj,
                    method = "gini",
                    expression_values = "normalized",
                    cluster_column = "new_celltype",
                    min_genes = 20,
                    min_expr_gini_score = 0.5,
                    min_det_gini_score = 0.5
                )
                print(head(gini_markers))
        
                # signature genes
                sign_markers = unique(
                    gini_markers[comb_rank <= 1000,genes]
                )
        
                # average expression per cell type
                average_cell_type_expr = Giotto:::create_average_DT(
                    gobject = sc_obj,
                    meta_data_name = "new_celltype",
                    expression_values = "normalized"
                )
        
                average_cell_type_expr = average_cell_type_expr[sign_markers, ]
        
                # ====== 官方示例里有、你原来缺的关键一步 ======
                colnames(average_cell_type_expr) =
                    gsub("cluster_", "", colnames(average_cell_type_expr))
    """)

        
    # Load spatial data into R
    print("Loading spatial data into R")
    scexp_st = anndata2ri._py2r.py2rpy_anndata(adata_spatial_copy)
    robjects.r.assign("scexp_st", scexp_st)
    robjects.r("""  
                counts_st = assay(scexp_st,"X") 
                meta_st = as.data.table(scexp_st@colData)
                meta_st$cell_ID = rownames(scexp_st@colData)
                gene_meta_st = as.data.table(rowData(scexp_st))
                gene_meta_st$feat_ID = rownames(rowData(scexp_st))
                coords = as.data.frame(scexp_st@int_colData@listData$reducedDims$spatial)
                colnames(coords) = c("x", "y")
                rownames(coords) = colnames(counts_st)
                giotto_obj_spatial = createGiottoObject(
                  raw_exprs     = counts_st,
                  spatial_locs  = coords,
                  cell_metadata = meta_st,
                  gene_metadata = gene_meta_st,
                  instructions  = instrs
                )
                giotto_obj_spatial = filterGiotto(gobject = giotto_obj_spatial, expression_threshold = 1, gene_det_in_min_cells = 5,
                      min_det_genes_per_cell = 5, expression_values = c('raw'), verbose = T)
                giotto_obj_spatial = normalizeGiotto(gobject = giotto_obj_spatial, scalefactor = 6000, verbose = T)
                giotto_obj_spatial <- addStatistics(gobject = giotto_obj_spatial)

                """)

    # Run deconvolution 
    print("Running deconvolution")
    robjects.r.assign("n_cell", n_cell)
    robjects.r.assign("cluster_col", cluster_key)
    robjects.r("giotto_obj_spatial = runDWLSDeconv(giotto_obj_spatial, cluster_column=cluster_col, sign_matrix=average_cell_type_expr, n_cell=n_cell)")

    # Save results 
    print("Saving results")
    robjects.r.assign("results_path", results_path)
    robjects.r("""
                dir.create(results_path)
                results_dt= giotto_obj_spatial@spatial_enrichment$DWLS
                write.csv(results_dt, paste0(results_path, '/proportions.csv'),row.names = TRUE)
                """)