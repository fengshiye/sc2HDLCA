def RCTD(
    adata_spatial, 
    adata_ref, 
    labels_key,
    r_lib_path=None, 
    results_path=None, 
): 
    from rpy2 import robjects
    import anndata2ri
    import scanpy as sc
    ###
    results_path=f'{results_path}/rctd_results/'
    adata_spatial_copy = adata_spatial.copy()
    adata_ref_copy = adata_ref.copy()
    scexp_sc = anndata2ri._py2r.py2rpy_anndata(adata_ref_copy)
    robjects.r.assign("scexp_sc", scexp_sc)
    scexp_st = anndata2ri._py2r.py2rpy_anndata(adata_spatial_copy)
    robjects.r.assign("scexp_st", scexp_st)
    ###3
    if r_lib_path is not None:
        robjects.r.assign("r_path", r_lib_path)
        robjects.r(".libPaths(r_path)")
    ###
    robjects.r.assign("labels_key", labels_key)
    robjects.r.assign("results_path", results_path)
    robjects.r("""
                suppressPackageStartupMessages({
                    library(spacexr)
                    library(Matrix)
                    library(data.table)
                    library(SingleCellExperiment) 
                })
                meta_sc=as.data.table(colData(scexp_sc))
                cell_types = as.vector(meta_sc[[labels_key]])
                names(cell_types) <- rownames(colData(scexp_sc))
                cell_types = factor(cell_types)
                sc_reference=Reference(
                    counts=assay(scexp_sc,'X'),
                    cell_types=cell_types
                )
                ###
                counts_st=assay(scexp_st,'X')
                location=as.data.frame(scexp_st@int_colData@listData$reducedDims$spatial)
                colnames(location)=c('x','y')
                rownames(location)=colnames(counts_st)
                st_data=SpatialRNA(
                    counts=counts_st,
                    coords=location,
                    require_int=FALSE
                )
                myRCTD <- create.RCTD(
                    spatialRNA = st_data,
                    reference = sc_reference,
                    max_cores = 1
                )
    
                myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
                ###
                dir.create(results_path, recursive = TRUE, showWarnings = FALSE)
                saveRDS(
                    myRCTD@results,
                    file = file.path(results_path, "RCTD_results.rds")
                )
            """)