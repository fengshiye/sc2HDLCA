from rpy2 import robjects
from rpy2.robjects import pandas2ri
import os
import scanpy as sc
from scvi.model import SCANVI
import torch
from rpy2.robjects import pandas2ri, conversion
from sc2hdlca.algorithms.sccafd import run_SCCAF_D  
from sc2hdlca.algorithms.popv import run_popv  
from sc2hdlca.algorithms.rctd import RCTD
from sc2hdlca.algorithms.scarches import train_scarches_scvi_scanvi_model
from sc2hdlca.algorithms.scarches import transfer_scarches_labels
from sc2hdlca.algorithms.scvi_scanvi import train_scvi_scanvi_model
from sc2hdlca.algorithms.scvi_scanvi import transfer_scANVI_labels
from sc2hdlca.algorithms.spatialdwls import spatialdwls
from sc2hdlca.algorithms.uniport_atac import run_uniPort
from sc2hdlca.algorithms.cellmarkeraccordion import run_cellmarkeraccordion
####
from importlib import resources

def get_example_adata():
    with resources.path("sc2hdlca.data", "HDLCA_reference.h5ad") as fpath:
        adata = sc.read_h5ad(fpath)
    return adata
###
def run_scHDLCA(
    adata_spatial=None,
    adata_query=None,
    adata_atac=None,
    labels_key='CellType_HDLCA',
    cluster_key=None,
    n_cell=50,
    deconvolution_type=None,
    at_sc=False,
    run_label_transfer=True,
    run_projection=True,
    run_annotation=True,
    results_path="./scHDLCA_results",
    python_home="",
    R_home="",
    bulk_path=None,
    query_batch_key=None,
    ref_batch_key=None
):
    """
    Main pipeline for scHDLCA framework.

    This function integrates spatial transcriptomics, single-cell RNA-seq,
    and optionally scATAC-seq data to perform:
    
    - Cell type deconvolution
    - Label transfer
    - Projection
    - Cell type annotation

    Parameters
    ----------
    adata_spatial : anndata.AnnData, optional
        Spatial transcriptomics dataset.

    adata_query : anndata.AnnData, optional
        Single-cell RNA-seq dataset used as query.

    adata_atac : anndata.AnnData, optional
        Single-cell ATAC-seq dataset.

    labels_key : str, optional
        Column name in reference AnnData (`adata_query.obs`) indicating cell type labels.

    cluster_key : str, optional
        Column name in `adata_spatial.obs` representing cluster assignments.

    n_cell : int, default=50
        Estimated number of cells per spatial spot.

    deconvolution_type : {"bulk", "spatial", None}, optional
        Type of cell type deconvolution to perform:
        - "bulk": bulk RNA-based deconvolution (via R)
        - "spatial": spatial-based deconvolution
        - None: skip deconvolution

    at_sc : bool, default=False
        Whether to integrate scATAC-seq data with reference data.

    run_label_transfer : bool, default=True
        Whether to perform label transfer from reference to query data.

    run_projection : bool, default=True
        Whether to project query data into reference space.

    run_annotation : bool, default=True
        Whether to perform cell type annotation.

    results_path : str, default="./scHDLCA_results"
        Directory to save output results.

    python_home : str, optional
        Path to Python executable used for bulk deconvolution (if required).

    R_home : str, optional
        Path to Rscript executable used for bulk deconvolution.

    bulk_path : str, optional
        Path to bulk RNA-seq data (required if `deconvolution_type="bulk"`).

    query_batch_key : str, optional
        Batch key in `adata_query.obs`.

    ref_batch_key : str, optional
        Batch key in reference dataset (`adata_query.obs` or others).

    Returns
    -------
    None
        Results are saved to `results_path`.

    Notes
    -----
    - This function may call external R scripts via rpy2 or subprocess.
    - Ensure that both Python and R environments are correctly configured.
    - For bulk deconvolution, `bulk_path`, `python_home`, and `R_home` must be provided.

    Examples
    --------
    >>> run_scHDLCA(
    ...     adata_spatial=adata_sp,
    ...     adata_query=adata_sc,
    ...     labels_key="celltype",
    ...     cluster_key="leiden",
    ...     deconvolution_type="bulk",
    ...     bulk_path="bulk.rds"
    ... )
    """

    # # 确保 R 脚本路径正确
    # r_dir = Path(__file__).parent.parent / "R/SCCAF-D"
    
    # # 设置 R 工作目录
    # robjects.r(f'setwd("{r_dir}")')

    # # 调用 R 脚本
    # robjects.r('source("SCCAF_D.R")')
    # robjects.r('source("benchmark1.R")')
    # robjects.r('source("deconvolution1.R")')
    # robjects.r('source("Frame.R")')
    # robjects.r('source("DWLS.R")')
    # robjects.r('source("expr.R")')
    ###
    print("sc2HDLCA is running...")
    os.makedirs(results_path, exist_ok=True)

    # ===== 直接加载 reference（固定）=====
    print("Loading reference dataset...")
    adata_ref = get_example_adata()
    print(adata_ref.shape)

    # ---------- 传对象给 R ----------
    if labels_key is not None:
        robjects.globalenv["labels_key"] = labels_key

    if cluster_key is not None:
        robjects.globalenv["cluster_key"] = cluster_key

    if n_cell is not None:
        robjects.globalenv["n_cell"] = n_cell

    if results_path is not None:
        robjects.globalenv["results_path"] = results_path

    if run_annotation:
        # ---------- Marker-based ----------
        print("Run Marker-based method (cellmarkeraccordion)")
        run_cellmarkeraccordion(
            adata_query,
            cluster_key,
            results_path=os.path.join(results_path, "cellmarkeraccordion")
        )

        # ---------- Reference-based ----------
        print("Run reference-based methods")
        run_popv(
            adata_query,
            adata_ref,
            labels_key,
            query_batch_key,
            ref_batch_key,
            results_path
        )

    # ---------- Deconvolution ----------
    if deconvolution_type is not None:
        if deconvolution_type == "bulk":
            print("Run Cell type deconvolution [bulk]")
            if bulk_path is None:
                raise ValueError("bulk_path must be specified when deconvolution_type='bulk'")
            run_SCCAF_D(
                bulk=bulk_path,
                python_home=python_home,
                results_path=os.path.join(results_path, "bulk_deconvolution"),
                R_home=R_home
            )
        elif deconvolution_type == "spatial":
            print("Run Cell type deconvolution [spatial]")
            print("Run RCTD!")
            RCTD(
                adata_spatial,
                adata_ref,
                labels_key,
                results_path=os.path.join(results_path, "rctd_results")
            )
            print("Run SpatialDWLS!")
            spatialdwls(
                adata_spatial,
                adata_ref,
                labels_key=labels_key,
                cluster_key=cluster_key,
                n_cell=n_cell,
                python_path=python_home,
                results_path=os.path.join(results_path, "spatialdwls_results")
            )
        else:
            raise ValueError("deconvolution_type must be 'bulk' or 'spatial'")

    # ---------- Label-aware integration ----------
    if run_label_transfer:
        print("Run label-aware transfer (scVI / scANVI)")
        out = os.path.join(results_path, "scvi_scanvi_results")
        os.makedirs(out, exist_ok=True)
        ###
        ad, model, scanvi_model = train_scvi_scanvi_model(
            adata_ref,
            seed=0,
            scvi_latent_key="X_scVI",
            batch_key="sample_ID",
            layer="counts",
            n_layers=1,
            n_hidden=128,
            n_latent=10,
            gene_likelihood="nb",
            labels_key=labels_key,
            unlabeled_category="Unlabeled",
            scanvi_latent_key="X_scANVI",
            scanvi_predictions_key="C_scANVI",
            save_model_path=results_path
        )
        ####
        adata_query, scanvi_query=transfer_scANVI_labels(
            adata_query,
            max_epochs=5,
            prediction_key="predicted_labels",
            latent_key="X_scANVI_query",
            results_path=results_path
        )
        
        adata_query.write(f'{out}/adata_query.h5ad')

    # ---------- Projection ----------
    if run_projection:
        print("Run projection (scArches)")
        out = os.path.join(results_path, "scarches_results")
        os.makedirs(out, exist_ok=True)
        ad, model, scanvi_model = train_scarches_scvi_scanvi_model(
            adata_ref,
            seed=0,
            scvi_latent_key="X_scVI_scArches",
            batch_key="sample_ID",
            layer="counts",
            n_layers=2,
            n_hidden=128,
            n_latent=10,
            gene_likelihood="nb",
            labels_key=labels_key,
            unlabeled_category="unlabel",
            scanvi_latent_key="X_scANVI_scArches",
            scanvi_predictions_key="C_scANVI_scArches",
            save_model_path=results_path
        )

        adata_query, scanvi_query=transfer_scarches_labels(
            adata_query,
            max_epochs=5,
            prediction_key="predicted_labels_scarches",
            latent_key="X_scarches_query",
            results_path=results_path
        )
        adata_query.write(f'{out}/adata_query.h5ad')

    # ---------- Optional ATAC sc integration ----------
    if at_sc:
        print('Run uniPort!')
        run_uniPort(
            adata_atac=adata_atac,
            adata_ref=adata_ref,
            output_path=os.path.join(results_path, "uniPort"),
            prefix="AT_sc_result"
        )

    print("scHDLCA pipeline finished.")