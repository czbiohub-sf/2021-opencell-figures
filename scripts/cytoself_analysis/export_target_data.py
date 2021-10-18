
'''
This script calculates the target labels dataframe and the anndata object for the target vectors
that are required by the umap viewer app

It also documents the workflow required to reproduce the raw target vectors, 
the leiden clusters, and the target UMAP coordinates
'''

from dataclasses import dataclass
import os
import scanpy as sc
import numpy as np
from . import model_results, clustering_workflows, ground_truth_labels


@dataclass
class Config:
    root_dirpath: str
    do_leiden_clustering: bool


def get_config(env):
    '''
    env: 'local' or 'ess'
    '''
    if env == 'local':
        config = Config(
            root_dirpath='/Users/keith.cheveralls/clustering-results', 
            do_leiden_clustering=True
        )

    # leiden clustering does not work on ESS (because the leidenalg package doesn't work)
    if env == 'ess':
        config = Config(
            root_dirpath='/ML_group/opencell-microscopy/clustering-results',
            do_leiden_clustering=False
        )
    return config


def main(env):

    config = get_config(env)

    # model is 'full', 'no-channel-splitting', or 'no-decoder'
    res = model_results.ModelResults.load_december_results(
        root_dirpath=config.root_dirpath,
        dataset='full', 
        model='no-decoder',
        rep=5
    )

    # concatenate the orphan test labels, test data, and VQ2 vectors 
    # res.concatenate_orphans()

    # append metadata fields and ground-truth labels
    res.test_labels = ground_truth_labels.merge_all(df=res.test_labels, data_dirpath='./data')

    if os.path.exists(res.get_target_labels_filepath()):
        print(
            'Error: cannot export the target labels CSV because one already exists in %s' 
            % res.exports_dirpath
        )
        return

    adata = res.export_adata(vq='vq2', kind='vectors', using='mean', pub_ready_only=True)
    cw = clustering_workflows.ClusteringWorkflow(adata=adata)
    
    # calculate neighbors
    cw.preprocess(
        do_log1p=False,
        do_scaling=False,
        n_top_genes=None,
        n_pcs=200,
    )
    cw.calculate_neighbors(n_neighbors=10, n_pcs=200, metric='euclidean')

    # reset the umap to the canonical parameters
    sc.tl.umap(
        cw.adata, 
        init_pos='spectral', 
        min_dist=0.0,
        random_state=42
    )

    if config.do_leiden_clustering:
        # a range of resolutions
        for resolution in [1, 2, 3, 5, 10, 20, 30]:
            cw.run_leiden(
                resolution, random_state=42, key_added='cluster_id_leiden_res=%0.1f' % resolution
            )
        # a range of random seeds at resolution=30
        for seed in range(10, 21):
            cw.run_leiden(
                resolution=30, 
                random_state=seed, 
                key_added='cluster_id_leiden_seed%s' % seed
            )

    # reset the umap to the canonical parameters
    sc.tl.umap(
        cw.adata, 
        init_pos='spectral', 
        min_dist=0.0,
        random_state=42
    )

    print('Exporting the target labels dataframe')
    umap_coords = cw.adata.obsm['X_umap']
    target_labels = cw.adata.obs.copy()
    target_labels['umap_0'] = umap_coords[:, 0]
    target_labels['umap_1'] = umap_coords[:, 1]
    target_labels.to_csv(res.get_target_labels_filepath())
    print('Target labels saved to %s' % res.get_target_labels_filepath())

    print('Exporting the adata object for the target VQ2 vectors')
    filepath = os.path.join(res.exports_dirpath, 'vq2-target-vectors-adata.h5ad')
    cw.adata.write_h5ad(filepath)


if __name__ == '__main__':
    main(env='local')
