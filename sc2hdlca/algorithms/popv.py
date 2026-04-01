def run_popv(
    adata_query, 
    adata_ref, 
    labels_key, 
    query_batch_key=None,
    ref_batch_key=None,
    results_path=None
):
    def _safe_predict(self, adata):
        from celltypist import train, models
        import os
        import pickle
    # 1. 转换 mask 为 numpy（核心修复！）
        mask = adata.obs["_predict_cells"] == "relabel"
        mask_np = mask.to_numpy()  # or mask.values
        
        # 2. 手动重建 AnnData（不依赖 view.copy()）
        sub = sc.AnnData(
            X=adata.X[mask_np].copy(),  # 用 numpy mask 索引
            obs=adata.obs[mask].copy(),
            var=adata.var.copy(),
            uns=adata.uns.copy(),
        )
        if "X_pca" in adata.obsm:
            sub.obsm["X_pca"] = adata.obsm["X_pca"][mask_np].copy()

        adata_ref1 = adata_ref.copy()  # 保留原始计数

        # 2. 归一化
        sc.pp.normalize_total(adata_ref1, target_sum=1e4)
        
        # 3. 对数据做 log1p
        sc.pp.log1p(adata_ref1)
        
        # 4. 可选：确保 var_names 唯一
        adata_ref1.var_names_make_unique()
        adata_ref1.obs_names_make_unique()
        # # 3. 复制必要的 uns
        # if '_save_path_trained_models' in adata.uns:
        #     sub.uns['_save_path_trained_models'] = adata.uns['_save_path_trained_models']
        # if 'ref_prediction_keys' in adata.uns:
        #     sub.uns['ref_prediction_keys'] = adata.uns['ref_prediction_keys']
        
        
        # 4. 执行预测
        # model_path = os.path.join(output_folder, "celltypist.pkl")
        # print(model_path)
        # os.makedirs(output_folder, exist_ok=True)
        # 如果模型不存在，先训练
        # if not os.path.exists(model_path):
        #     self.compute_integration(adata)
        
       # 1. 训练模型
       # 训练模型
        model = train(
            adata_ref1,         # 已 normalize+log1p
            labels=labels_key
        )
        
        #####
        predictions = celltypist.annotate(
            sub,
            model=model,
            **self.classifier_dict,
        )
        
        # 5. 写回结果
        out_col = "majority_voting" if "majority_voting" in predictions.predicted_labels.columns else "predicted_labels"
        adata.obs.loc[mask, self.result_key] = predictions.predicted_labels[out_col].values
        
        return self
    ####
    import popv
    popv.settings.n_jobs = 10
    from popv.algorithms._celltypist import CELLTYPIST#
    import popv.algorithms._celltypist as ct_module#
    import scanpy as sc
    import numpy as np
    import popv
    import celltypist
    from celltypist import models
    import os
    
    _original_predict = CELLTYPIST.predict#
    CELLTYPIST.predict = _safe_predict#
    output_folder = f'{results_path}/scHDLCA_popv_model'
    os.makedirs(output_folder, exist_ok=True)
    unknown_celltype_label = "unassigned"  # Label of unlabeled cells
    n_samples_per_label = 300 # Downsamples for some classifiers the dataset.
    adata = popv.preprocessing.Process_Query(
    adata_query,
    adata_ref,
    query_batch_key=query_batch_key,
    ref_labels_key=labels_key,
    ref_batch_key=ref_batch_key,
    unknown_celltype_label=unknown_celltype_label,
    save_path_trained_models=output_folder,
    cl_obo_folder="/home/ubuntu/data/fengshuo/hdd130/fs/fs/ST/data/popv/ontology" ,
    prediction_mode="retrain",
    n_samples_per_label=n_samples_per_label,
    hvg=4000,
 ).adata.copy()

    ###
    # 关键：彻底 materialize
    adata = adata.to_memory() #scanvi error
    #
    adata.var_names_make_unique()
    adata.obs_names_make_unique() 
    ##
    # 完全关闭 ontology graph 相关逻辑
    popv.annotation.ontology_vote_onclass = lambda *args, **kwargs: None
    popv.annotation.ontology_parent_onclass = lambda *args, **kwargs: None
    ###
    for col in [
        "popv_prediction",
        "popv_prediction_score",
        "popv_parent",
    ]:
        adata.obs[col] = np.nan
    ###
    popv.annotation.annotate_data(
        adata,methods=[
            'SCANVI_POPV',
            "CELLTYPIST",
            "Support_Vector",
            "XGboost",
            "Random_Forest",
            "KNN_SCANORAMA",
            "KNN_SCVI"
        ],
        save_path=f"{results_path}/popv_output",
    )