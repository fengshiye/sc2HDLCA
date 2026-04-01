import scvi
import torch
from tqdm import tqdm
import os
def train_scvi_scanvi_model(
    adata,  # AnnData 对象
    seed=0,
    scvi_latent_key="X_scVI",
    batch_key='sample_ID',
    continuous_covariate_keys=None,
    layer="counts",
    n_layers=1,
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
    max_epochs=None,  # 默认为None，如果设置为None，则使用early stopping
    save_model_path=None,  # 模型保存路径
    labels_key="celltype_combined_coarse",
    unlabeled_category="Unlabeled",
    scanvi_latent_key="X_scANVI_sample",
    scanvi_predictions_key="C_scANVI_sample",
):  
    # 设置随机种子
    scvi.settings.seed = seed

    # 获取高变
    adataHVG = adata[:,adata.var['highly_variable']].copy()
    
    # 设置 AnnData 对象用于 scVI 模型
    scvi.model.SCVI.setup_anndata(
        adataHVG,
        layer=layer,
        batch_key=batch_key,
        continuous_covariate_keys=continuous_covariate_keys
    )

    # 初始化并训练 scVI 模型
    model = scvi.model.SCVI(
        adataHVG,
        n_layers=n_layers,
        n_hidden=n_hidden,
        n_latent=n_latent,
        dispersion=dispersion,
        gene_likelihood=gene_likelihood,
        encode_covariates=encode_covariates,
        deeply_inject_covariates=deeply_inject_covariates,
        use_layer_norm=use_layer_norm,
        use_batch_norm=use_batch_norm,
    )
    model.train(
        check_val_every_n_epoch=check_val_every_n_epoch,
        early_stopping=early_stopping,
        early_stopping_patience=early_stopping_patience,
        early_stopping_monitor=early_stopping_monitor,
        batch_size=batch_size,
        train_size=train_size,
        max_epochs=max_epochs,
    )
    
    # # 绘制训练和验证的重建损失曲线
    # plt.plot(model.history['reconstruction_loss_train']['reconstruction_loss_train'], label='train')
    # plt.plot(model.history['reconstruction_loss_validation']['reconstruction_loss_validation'], label='validation')
    # plt.legend()
    # plt.show()

    # 保存模型（如果提供了保存路径）
    # if save_model_path is not None:
    #     model.save(save_model_path, overwrite=True)

    # 使用 scVI 模型创建 scANVI 模型
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        model,
        adata=adataHVG,
        labels_key=labels_key,
        unlabeled_category=unlabeled_category,
    )
    scanvi_model.train()
    ###
    if save_model_path is not None:
        model.save(os.path.join(save_model_path, "scvi_scanvi_results",'scvi'), overwrite=True)
        scanvi_model.save(os.path.join(save_model_path, "scvi_scanvi_results",'scanvi'), overwrite=True)
        print('save model...')
    # # 绘制 scANVI 模型的训练损失曲线
    # plt.plot(scanvi_model.history['reconstruction_loss_train']['reconstruction_loss_train'], label='train')
    # plt.legend()
    # plt.show()

    # 保存 scANVI 模型的潜在表示
    adata.obsm[scanvi_latent_key] = scanvi_model.get_latent_representation(adataHVG)

    # 预测细胞类型并保存到 AnnData 对象
    adata.obs[scanvi_predictions_key] = scanvi_model.predict(adataHVG)

    return adata, model, scanvi_model
####
def transfer_scANVI_labels(
    adata_query,
    max_epochs=5,
    prediction_key="predicted_labels_scanvi",
    latent_key="X_scANVI_query",
    results_path=None
    
):
   ref_path=os.path.join(results_path, "scvi_scanvi_results",'scanvi')
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
