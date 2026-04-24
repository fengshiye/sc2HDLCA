# 🧬 sc2HDLCA

**sc2HDLCA** is a unified toolkit designed to facilitate reference-based
understanding of lung single-cell biology.

It leverages the **HDLCA (Human Developmental Lung Cell Atlas)** as a
reference to enable comprehensive analysis across multiple data
modalities.

------------------------------------------------------------------------

## 🚀 Overview

To facilitate reference-based analysis of lung single-cell data,
**sc2HDLCA** provides a flexible framework for:

-   🏷️ **Cell type annotation** of new dataset
-   🔬 **Cell type deconvolution across modalities**
-   🔗 **Integration and analysis of multi-omics data**, including:
    -   Single-cell RNA-seq
    -   Spatial transcriptomics
    -   Bulk RNA-seq

------------------------------------------------------------------------

## ✨ Key Features

-   **Reference-based annotation**
    Annotate new datasets using HDLCA as a high-quality reference atlas

-   **Cross-modality deconvolution**
    Perform cell type deconvolution for:

    -   Bulk RNA-seq
    -   Spatial transcriptomics

-   **Multi-omics integration**
    Support integration of:

    -   scRNA-seq
    -   scATAC-seq (optional)

-   **Modular design**
    Flexible pipeline where each component can be run independently

------------------------------------------------------------------------

## 📦 Installation sc2HDLCA
Firstly, download the environment file and reference data [here](https://your-link-here.com).<br></br>
Secondly, install the environment.<br></br>
Example:
``` bash
  mkdir sc2hdlca_env
  tar -xzf sc2hdlca_env.tar.gz -C sc2hdlca_env
  cd sc2hdlca_env
  ./bin/conda-unpack
```

------------------------------------------------------------------------

## 🧪 Usage

``` python
from sc2hdlca.algorithms.scHDLCA_core import run_scHDLCA

run_scHDLCA(
    adata_spatial=adata_sp,
    adata_query=adata_sc,
    labels_key="celltype",
    cluster_key="leiden",
    deconvolution_type="bulk",
    bulk_path="bulk_data.rds",
    results_path="./results"
)
```
## Parameters

### Cell Annotation
- **adata_query**: Single-cell RNA-seq dataset used as query  
- **run_label_transfer**: Whether to perform label transfer from reference to query data  

---

### Cell Annotation
- **adata_query**: Single-cell RNA-seq dataset used as query  
- **run_annotation**: Whether to perform cell type annotation  

---

### Bulk Data Deconvolution
- **bulk_path**: Path to bulk RNA-seq data (required if `deconvolution_type="bulk"`)  
- **R_home**: Path to Rscript executable used for bulk deconvolution  
- **python_home**: Path to Python executable used for bulk deconvolution (if required)  

---

### Spatial Data Deconvolution
- **adata_spatial**: Spatial transcriptomics dataset  
- **cluster_key**: Column name in `adata_spatial.obs` representing cluster assignments  
- **deconvolution_type**: Type of deconvolution to perform:
  - `"bulk"`: bulk RNA-based deconvolution (via R)  
  - `"spatial"`: spatial-based deconvolution  
  - `None`: skip deconvolution  

---

### scATAC Integration
- **adata_atac**: Single-cell ATAC-seq dataset  
- **at_sc**: Whether to integrate scATAC-seq data with reference data  

---

### Projection
- **adata_query**: Single-cell RNA-seq dataset used as query  
- **run_projection**: Whether to project query data into reference space  
------------------------------------------------------------------------

## 📁 Input Data

-   `adata_query`: single-cell RNA-seq
-   `adata_spatial`: spatial transcriptomics
-   `bulk_path`: bulk RNA-seq
-   `adata_atac`: single-cell ATAC-seq

------------------------------------------------------------------------

## 📚 Citation

    [Your paper / preprint here]

------------------------------------------------------------------------

## 📬 Contact

For questions or collaborations, please open an issue.
