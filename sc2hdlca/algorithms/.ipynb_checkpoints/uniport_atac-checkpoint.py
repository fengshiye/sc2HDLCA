# %%
import os
import sys
import uniport as up
import numpy as np
import pandas as pd
import torch
import scanpy as sc
import episcanpy as epi
from sklearn.preprocessing import MinMaxScaler
import scipy.sparse
from nvitop import Device

# 自动选择空闲 GPU
devices = Device.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)


# %%
def run_uniPort(
    adata_atac=None,
    adata_ref=None,
    # input_path,
    output_path=None,
    # rna_file="rna.h5ad",
    # atac_file="atac.h5ad",
    prefix="result"
):
    # Make Dir
    os.makedirs(output_path, exist_ok=True)

    print("Loading data...")
    # rna = sc.read_h5ad(os.path.join(input_path, rna_file))
    # atac = sc.read_h5ad(os.path.join(input_path, atac_file))
    atac=adata_atac
    rna=adata_ref
    ####
    import scipy.sparse as sp

    if sp.issparse(rna.X):
        rna.X = rna.X.toarray()
    
    if sp.issparse(atac.X):
        atac.X = atac.X.toarray()
    # 设置 domain
    atac.obs["domain_id"] = 0
    atac.obs["domain_id"] = atac.obs["domain_id"].astype("category")
    atac.obs["source"] = "ATAC"

    rna.obs["domain_id"] = 1
    rna.obs["domain_id"] = rna.obs["domain_id"].astype("category")
    rna.obs["source"] = "RNA"

    # ================= RNA preprocess =================
    print("Preprocess RNA...")
    sc.pp.normalize_total(rna)
    sc.pp.log1p(rna)
    sc.pp.highly_variable_genes(rna, n_top_genes=3000, subset=True)
    up.batch_scale(rna, chunk_size=20001)

    # ================= ATAC preprocess =================
    print("Preprocess ATAC...")
    atac.X[atac.X > 1] = 1
    epi.pp.select_var_feature(atac, nb_features=2000, show=False, copy=False)
    sc.pp.normalize_total(atac)
    up.batch_scale(atac, chunk_size=20001)

    # ================= Run uniPort =================
    print("Running uniPort...")
    comb_tab = up.Run(
        adatas=[atac, rna],
        mode="v",
        lr=0.0001,
        iteration=10000,
        num_workers=False
    )

    # ================= Cluster =================
    print("Clustering...")
    sc.pp.neighbors(comb_tab, use_rep="latent")
    sc.tl.umap(comb_tab, min_dist=0.1)
    sc.tl.leiden(comb_tab)

    # ================= Save latent =================
    latent = pd.DataFrame(
        comb_tab.obsm["latent"],
        index=comb_tab.obs_names
    )
    latent.to_csv(
        os.path.join(output_path, f"{prefix}-uniPort-multi-latent.csv")
    )

    # ================= Save UMAP =================
    umap = pd.DataFrame(
        comb_tab.obsm["X_umap"],
        columns=["UMAP1", "UMAP2"],
        index=comb_tab.obs_names
    )
    umap.insert(2, "cluster", comb_tab.obs["leiden"].values)

    umap.to_csv(
        os.path.join(output_path, f"{prefix}-uniPort-multi-umap.csv")
    )

    # ================= GPU memory =================
    pid = os.getpid()
    gpu_memory = pd.Series(dtype="str")

    devices = Device.all()
    for device in devices:
        processes = device.processes()
        if pid in processes:
            p = processes[pid]
            gpu_memory[f"device {device.index}"] = p.gpu_memory_human()

    if len(gpu_memory):
        gpu_memory.to_csv(
            os.path.join(output_path, f"{prefix}-uniPort-gpu_memory.csv"),
            header=["gpu_memory"]
        )

    print("uniPort finished!")


# %%
# 支持命令行运行
# if __name__ == "__main__":
#     input_path = sys.argv[1]
#     output_path = sys.argv[2]

#     uniPort_module(
#         input_path=input_path,
#         output_path=output_path
#     )