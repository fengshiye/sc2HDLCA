def run_cellmarkeraccordion(
    adata_query,
    cluster_key,
    results_path
):
    from rpy2 import robjects
    import anndata2ri
    import scanpy as sc
    adata_query_copy = adata_query.copy()
    scexp_sc = anndata2ri._py2r.py2rpy_anndata(adata_query_copy)
    robjects.r.assign("scexp_sc", scexp_sc)
    robjects.r.assign("cluster_key", cluster_key)
    robjects.r.assign("results_path", results_path)
    robjects.r("""
            suppressPackageStartupMessages({
                library(cellmarkeraccordion)
                library(Seurat)
                library(data.table)
                library(ggplot2)
            })
            ###
            meta_sc=as.data.table(colData(scexp_sc))
            counts_sc=assay(scexp_sc,'X')
            ###
            # Create Seurat Object
            data <- CreateSeuratObject(counts = counts_sc,meta.data=meta_sc,min.cells = 3, min.features = 200)
            ##
            HDLCA_markers=read.csv('/home/ubuntu/data/fengshuo/project/devlung_2025-10-30/HDLCA_marker_gene2celltype.csv')#read cell marker table
            retinal_data <- accordion_custom(
                data,
                annotation_resolution = cluster_key,
                marker_table = HDLCA_markers,
                category_column= "terms",
                marker_column ="Symbol",
                min_n_marker = 2,
                plot = FALSE,
                annotation_name = "cell_type_HDLCA"
            )
            dir.create(results_path)
            saveRDS(
                retinal_data,
                file = file.path(results_path, "cellmarkeraccordion_results.rds")
            )
        """)