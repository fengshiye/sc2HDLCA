import scarches as sca
import scvi
import torch
import os
def train_scarches_scvi_scanvi_model(
    adata,  # AnnData 对象
    seed=0,
    scvi_latent_key="X_scVI_scArches",
    batch_key='sample_ID',
    continuous_covariate_keys=None,
    layer="counts",
    n_layers=2,
    n_hidden=128,
    n_latent=10,
    gene_likelihood="nb",
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
    dispersion="gene-batch",
    check_val_every_n_epoch=1,
    early_stopping=True,
    early_stopping_patience=20,
    early_stopping_monitor="elbo_validation",
    batch_size=740,
    train_size=0.75,
    max_epochs=None,  # 默认为None,如果设置为None,则使用early stopping
    save_model_path=None,  # 模型保存路径
    labels_key="celltype_final",
    unlabeled_category="unlabel",
    scanvi_latent_key="X_scANVI_scArches",
    scanvi_predictions_key="C_scANVI_scArches",
):
    """
    使用 scArches 训练 SCVI 和 SCANVI 模型用于参考数据
    
    Parameters:
    -----------
    adata : AnnData
        包含 'highly_variable' 列的 AnnData 对象
    其他参数与 SCANVI 函数相同
    
    Returns:
    --------
    adata : AnnData
        更新后的 AnnData 对象,包含潜在表示和预测结果
    vae : sca.models.SCVI
        训练好的 scVI 模型
    scanvae : sca.models.SCANVI
        训练好的 scANVI 模型
    """


    # 获取高变基因并准备数据
    adataHVG = adata[:, adata.var['highly_variable']].copy()
    adataHVG.X = adataHVG.layers[layer].copy()
    
    # 移除稀疏性(如果需要)
    # from scarches.dataset.trvae.data_handling import remove_sparsity
    # adataHVG = remove_sparsity(adataHVG)
    
    print(f"Training on {adataHVG.n_obs} cells and {adataHVG.n_vars} genes")
    
    # 设置 AnnData 对象用于 scArches SCVI 模型
    sca.models.SCVI.setup_anndata(
        adataHVG,
        batch_key=batch_key,
        labels_key=labels_key,
        layer=None,  # 因为已经复制到X了
        continuous_covariate_keys=continuous_covariate_keys
    )
    
    # 初始化并训练 scArches scVI 模型
    vae = sca.models.SCVI(
        adataHVG,
        n_layers=n_layers,
        n_hidden=n_hidden,
        n_latent=n_latent,
        gene_likelihood=gene_likelihood,
        encode_covariates=encode_covariates,
        deeply_inject_covariates=deeply_inject_covariates,
        use_layer_norm=use_layer_norm,
        use_batch_norm=use_batch_norm,
    )
    
    vae.train(
        max_epochs=max_epochs,
        check_val_every_n_epoch=check_val_every_n_epoch,
        early_stopping=early_stopping,
        early_stopping_patience=early_stopping_patience,
        early_stopping_monitor=early_stopping_monitor,
        batch_size=batch_size,
        train_size=train_size,
    )
    
    # # 绘制 scVI 训练和验证的重建损失曲线
    # print("Plotting scVI training history...")
    # import matplotlib.pyplot as plt

    
    # plt.figure(figsize=(10,5))
    
    # # 训练损失
    # plt.plot(vae.history["elbo_train"]["loss"], label="train", linewidth=2)
    
    # # 验证损失（如果有）
    # if "elbo_validation" in vae.history:
    #     plt.plot(vae.history["elbo_validation"]["loss"], label="validation", linewidth=2)

    # plt.xlabel('Epoch')
    # plt.ylabel('Reconstruction Loss')
    # plt.title('scVI Training History')
    # plt.legend()
    # plt.grid(True, alpha=0.3)
    # plt.tight_layout()
    # plt.show()
    
    # 保存 scVI 模型(如果提供了保存路径)
    if save_model_path is not None:
        scvi_model_path = save_model_path.replace('.pt', '_scvi.pt') if '.pt' in save_model_path else f"{save_model_path}_scvi"
        vae.save(scvi_model_path, overwrite=True)
        print(f"scVI model saved to: {scvi_model_path}")
    
    # 保存 scVI 潜在表示
    adata.obsm[scvi_latent_key] = vae.get_latent_representation(adataHVG)
    
    # 使用 scVI 模型创建 scANVI 模型
    print("\nInitializing scANVI model from scVI...")
    scanvae = sca.models.SCANVI.from_scvi_model(
        vae,
        unlabeled_category=unlabeled_category,
        adata=adataHVG,
        labels_key=labels_key,
    )
    
    print(f"Labeled Indices: {len(scanvae._labeled_indices)}")
    print(f"Unlabeled Indices: {len(scanvae._unlabeled_indices)}")
    
    # 训练 scANVI 模型
    scanvae.train(
        max_epochs=max_epochs,
        check_val_every_n_epoch=check_val_every_n_epoch,
        early_stopping=early_stopping,
        early_stopping_patience=early_stopping_patience,
        batch_size=batch_size,
        train_size=train_size,
    )
    
   #  # 绘制 scANVI 训练损失曲线
   #  print("Plotting scANVI training history...")
   #  plt.figure(figsize=(10, 5))
   # # 训练损失
   #  plt.plot(scanvae.history["elbo_train"]["loss"], label="train", linewidth=2)
    
   #  # 验证损失（如果有验证集）
   #  if "elbo_validation" in scanvae.history:
   #      plt.plot(scanvae.history["elbo_validation"]["loss"], label="validation", linewidth=2)

   #  plt.xlabel('Epoch')
   #  plt.ylabel('Reconstruction Loss')
   #  plt.title('scANVI Training History')
   #  plt.legend()
   #  plt.grid(True, alpha=0.3)
   #  plt.tight_layout()
   #  plt.show()
    
    # 保存 scANVI 模型(如果提供了保存路径)
    if save_model_path is not None:
        # scanvi_model_path = save_model_path.replace('.pt', '_scanvi.pt') if '.pt' in save_model_path else f"{save_model_path}_scanvi"
        # scanvae.save(scanvi_model_path, overwrite=True, save_anndata=True)
        vae.save(os.path.join(save_model_path, "scarches_results",'scvi'), overwrite=True)
        scanvae.save(os.path.join(save_model_path, "scarches_results",'scanvi'), overwrite=True)
        print(f"scANVI model saved ...")
    
    # 保存 scANVI 潜在表示
    adata.obsm[scanvi_latent_key] = scanvae.get_latent_representation(adataHVG)
    
    # 预测细胞类型并保存到 AnnData 对象
    adata.obs[scanvi_predictions_key] = scanvae.predict(adataHVG)
    
    print("\nTraining completed successfully!")
    print(f"scVI latent representation saved to: adata.obsm['{scvi_latent_key}']")
    print(f"scANVI latent representation saved to: adata.obsm['{scanvi_latent_key}']")
    print(f"Cell type predictions saved to: adata.obs['{scanvi_predictions_key}']")
    
    return adata, vae, scanvae

def transfer_scarches_labels(
    adata_query,
    max_epochs=5,
    prediction_key="predicted_labels_scarches",
    latent_key="X_scarches_query",
    results_path=None
):

   ref_path=os.path.join(results_path, "scarches_results",'scanvi')
    # 1. 对齐
   scvi.model.SCANVI.prepare_query_anndata(adata_query, ref_path)
    # 2. 加载
   scanvi_query=scvi.model.SCANVI.load_query_data(
        adata_query,
        ref_path,
        freeze_dropout = True,
            # inplace_subset_query_vars = True,
    )

    # 3. 微调
   scanvi_query.train(max_epochs=max_epochs)

    # 4. 预测
   adata_query.obs[prediction_key] = scanvi_query.predict()

    # 5. latent
   adata_query.obsm[latent_key] = scanvi_query.get_latent_representation()

   return adata_query, scanvi_query