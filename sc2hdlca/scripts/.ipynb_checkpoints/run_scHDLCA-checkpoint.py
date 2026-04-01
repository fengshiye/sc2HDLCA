#!/usr/bin/env python
from pathlib import Path
import sys

# 临时加入包路径，如果没 pip install
sys.path.append(str(Path(__file__).resolve().parents[1]))

from sc2hdlca.algorithms.scHDLCA_core import run_scHDLCA

if __name__ == "__main__":
    # 这里可以填你想用的参数
    run_scHDLCA(
    adata_spatial=None,
    adata_query=None,
    adata_atac=None,
    labels_key='celltype_final',
    cluster_key='leiden',
    n_cell=50,
    deconvolution_type='bulk',
    at_sc=False,
    run_label_transfer=False,
    run_projection=False,
    results_path="/home/ubuntu/data/fengshuo/hdd130/fs/fs/project/devlung/",
    python_home='/home/ubuntu/data/fengshuo/.envs/SCCAFD/bin/python',
    R_home='/home/ubuntu/data/fengshuo/.envs/SCCAFD/bin/Rscript',
    bulk_path='/home/ubuntu/data/fengshuo/hdd130/fs/fs/project/devlung/pseudobulk_Baron_T.rds',
    query_batch_key='sample_ID',
    ref_labels_key='celltype_final',
    ref_batch_key='sample_ID',
    run_annotation=False
    
)