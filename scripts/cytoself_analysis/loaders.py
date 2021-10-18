import os
import re
import io
import sys
import glob
import enum
import json
import umap
import dask
import natsort
import datetime
import skimage
import sklearn
import numpy as np
import scipy as sp
import pandas as pd
import seaborn as sns
import anndata as ad

from . model_results import ModelResults

MODEL_RESULTS_DIR = '/ML_group/opencell-microscopy/clustering-results'
GOOD_FOVS_FILEPATH = '/ML_group/opencell-microscopy/2020-05-27_good-fovs-with-labels.csv'


def maybe_load_data(filepath, loader):
    try:
        data = loader(filepath)
    except FileNotFoundError:
        print('Warning: file %s not found' % filepath)
        data = None
    return data


def load_october_results(root_dirpath=None):
    '''
    These October 2020 results are from a model with a VQ1 layer connected to an FC layer 
    that 'classifies' or 'tethers together' image patches from the same cell_line_id
    (note that there is no VQ2 layer in this model)
    '''
    nickname = 'october-results'
    if root_dirpath is None:
        root_dirpath = MODEL_RESULTS_DIR

    root_dirpath = os.path.join(root_dirpath, '2020-october-results')
    results_dirpath = os.path.join(root_dirpath, 'efb4zp_1lyfc_nd', 'fc2_1000_pretrain')

    test_label_array = np.load(
        os.path.join(root_dirpath, 'valtest_data', 'data_label.npy'),
        allow_pickle=True
    )

    # test_data only exists on ESS, not locally
    test_data = maybe_load_data(
        os.path.join(root_dirpath, 'valtest_data', 'data_image.npy'),
        loader=lambda filepath: np.load(filepath, allow_pickle=True)
    )

    # construct a dataframe from the numpy array of test labels
    # (note that the order of the column names is empirically determined to match this array)
    test_labels = pd.DataFrame(
        data=test_label_array[:, :13], 
        columns=[
            'cell_line_id', 'target_name', 'unknown_id', 'filepath',
            'label_0', 'label_1', 'label_2', 'label_3', 'label_4', 'label_5',
            'plate_id', 'well_id', 'pml_id',
        ]
    )

    # load the per-image UMAP coordinates (pre-computed and cached by Hiro)
    umap_coords = np.load(
        os.path.join(results_dirpath, 'umap_data', 'vt_vq1_indhist_umap.npy')
    )

    # load the VQ1 indices
    test_vq1_ind = np.load(
        os.path.join(results_dirpath, 'embeddings', 'vt_vq1_ind.npy')
    )

    # load the VQ1 vectors themselves (only on ESS)
    test_vq1 = maybe_load_data(
        os.path.join(results_dirpath, 'embeddings', 'vt_vq1_vec.npy'),
        loader=np.load
    )

    # append the UMAP coordinates to the patch labels dataframe
    test_labels['umap_0'] = umap_coords[:, 0]
    test_labels['umap_1'] = umap_coords[:, 1]

    # where to save exported results for the umap viewer
    exports_dirpath = os.path.join(results_dirpath, 'umap-viewer-data')

    return ModelResults(
        nickname=nickname,
        num_features=256,
        test_data=test_data,
        test_labels=test_labels,
        test_vq1=test_vq1,
        test_vq2=None,
        test_vq1_ind=test_vq1_ind,
        test_vq2_ind=None,
        exports_dirpath=exports_dirpath
    )


def load_december_results(dataset='full', model='full', rep=1, root_dirpath=None):
    '''
    These December 2020 results are from a model with VQ1 and VQ2 layers,
    both with an FC layer, and with multiple (four) 'channels' in the VQ2 layer 

    Two different results can be loaded: 
        1) if dataset == 'old', the 'initial' results from this model architecture,
            using our original training data defined in '2020-05-27_good-fovs.csv'
        2) if dataset == 'full', the latest results from this same architecture,
            but using the new training data defined in '2020-12-20_good-fovs.csv', 
            and also with crops centered on the nuclei

    Different models can be loaded:
        1) model == 'full' is the full model with all features (FC layers and channel splitting)
        2) model == 'no-channel-splitting' is the model without channel splitting
        3) model == 'no-decoder' is the model without the decoder
    '''

    if root_dirpath is None:
        root_dirpath = MODEL_RESULTS_DIR

    if dataset not in ['old', 'full']:
        raise ValueError

    if dataset == 'old':
        subdirname = 'olddata'
        appendix = ''
        test_label_columns = [
            'cell_line_id', 'target_name', 'unknown_id', 'filepath',
            'label_0', 'label_1', 'label_2', 'label_3', 'label_4', 'label_5',
            'plate_id', 'well_id', 'pml_id',
        ]

    elif dataset == 'full':
        subdirname = 'fulldata'
        appendix = '_nucenter'
        test_label_columns = [
            'cell_line_id', 'target_name', 'fov_id', 'filepath', 'plate_id', 'well_id', 'pml_id'
        ]

    nickname = 'december-results-%s' % dataset

    # the directory containing the results directory and the test data directory
    root_dirpath = os.path.join(root_dirpath, '2020-december-results', subdirname)

    # the name of the results directory 
    # note that for the 'full' dataset, we use the model trained with nucleus-centered crops
    results_dirname = '[gfp, nucdist]%s' % appendix

    # the results subdirectory itself
    if model == 'full':
        results_dirpath = os.path.join(
            root_dirpath,
            'efb0_chsplit_fc_nd2nd_trgt', 
            'z11[25, 25, 64]2048_z29[4, 4, 64]2048', 
            'fc2_1000_vq[1, 2]',
            results_dirname,
            'rep%d' % rep,
        )
        # this is the highest-scoring patch UMAP replicate from the random seeds that Hiro tried
        # (these seeds yield different UMAP 'replicates', which are separate 
        # from the model replicate indicated by the `rep%d` subdir name above)
        patch_umap_filename = 'urep6_vec2.npy'

    if model == 'no-channel-splitting':
        # only one rep is available for this model
        rep = 5
        results_dirpath = os.path.join(
            root_dirpath,
            'efb0zp_fc_nd2nd_trgt', 
            'z1[25, 25, 64]2048_z2[4, 4, 64]2048', 
            'fc2_1000_vq[1, 2]',
            results_dirname,
            'rep%d' % rep,
        )
        patch_umap_filename = 'urep5_vec2.npy'

    if model == 'no-decoder':
        # only one rep is available for this model
        rep = 1
        results_dirpath = os.path.join(
            root_dirpath,
            'encoder_only', 
            'z1[25, 25, 64]_z2[4, 4, 576]', 
            'fc2_1000_do[0.5, 0.5]',
            results_dirname,
            'rep%d' % rep,
        )
        patch_umap_filename = 'test_vqvec_umap2.npy'

    # the number of features in the codebook        
    num_features = 2048

    test_label_array = np.load(
        os.path.join(root_dirpath, 'test_data', 'test_label%s.npy' % appendix),
        allow_pickle=True
    )

    # test_data (the image patches themselves) only exists on ESS, not locally
    test_data = maybe_load_data(
        os.path.join(root_dirpath, 'test_data', 'test_data%s.npy' % appendix),
        loader=lambda filepath: np.load(filepath, allow_pickle=True)
    )

    # construct a dataframe from the numpy array of test labels
    # (note that the order of the column names is empirically determined to match this array)
    test_labels = pd.DataFrame(data=test_label_array[:, :13], columns=test_label_columns)

    # load the indices
    test_vq1_ind = maybe_load_data(
        os.path.join(results_dirpath, 'embeddings', 'test_vqind1.npy'),
        loader=np.load
    )
    test_vq2_ind = maybe_load_data(
        os.path.join(results_dirpath, 'embeddings', 'test_vqind2.npy'),
        loader=np.load
    )

    # load the vectors themselves (only on ESS)
    test_vq1 = maybe_load_data(
        os.path.join(results_dirpath, 'embeddings', 'test_vqvec1.npy'),
        loader=np.load
    )
    test_vq2 = maybe_load_data(
        os.path.join(results_dirpath, 'embeddings', 'test_vqvec2.npy'),
        loader=np.load
    )

    # load the 'orphan' vectors and test labels (these only exist for 'fulldata' replicates)
    orphan_vq1_ind = maybe_load_data(
            os.path.join(results_dirpath, 'embeddings', 'orphan_vqind1.npy'),
        loader=np.load
    )
    orphan_vq2_ind = maybe_load_data(
            os.path.join(results_dirpath, 'embeddings', 'orphan_vqind2.npy'),
        loader=np.load
    )
    orphan_vq2 = maybe_load_data(
            os.path.join(results_dirpath, 'embeddings', 'orphan_vqvec2.npy'),
        loader=np.load
    )
    orphan_label_array = np.load(
        os.path.join(root_dirpath, 'test_data', 'orphan_label%s.npy' % appendix),
        allow_pickle=True
    )
    orphan_data = maybe_load_data(
        os.path.join(root_dirpath, 'test_data', 'orphan_data%s.npy' % appendix),
        loader=lambda filepath: np.load(filepath, allow_pickle=True)
    )

    # load the per-patch UMAP coordinates from the VQ2 *vectors* (not histograms!) 
    # (pre-computed and cached by Hiro)
    umap_coords = np.load(
        os.path.join(results_dirpath, 'umap_data', patch_umap_filename)
    )

    # append the UMAP coordinates to the patch labels dataframe
    test_labels['umap_0'] = umap_coords[:, 0]
    test_labels['umap_1'] = umap_coords[:, 1]

    # where to save exported results for the umap viewer
    exports_dirpath = os.path.join(results_dirpath, 'umap-viewer-data')

    model_results = ModelResults(
        nickname=nickname,
        num_features=num_features,
        test_data=test_data,
        test_labels=test_labels,
        test_vq1=test_vq1,
        test_vq2=test_vq2,
        test_vq1_ind=test_vq1_ind,
        test_vq2_ind=test_vq2_ind,
        exports_dirpath=exports_dirpath
    )

    # addendum: load the pre-computed cached codebooks 
    try:
        model_results.codebook_vq1 = np.load(
            os.path.join(results_dirpath, 'embeddings', 'codebook_vq1.npy')
        )
        model_results.codebook_vq2 = np.load(
            os.path.join(results_dirpath, 'embeddings', 'codebook_vq2.npy')
        )
        # transpose to get (n_features, n_feature_dims)
        model_results.codebook_vq1 = model_results.codebook_vq1.transpose()
        model_results.codebook_vq2 = model_results.codebook_vq2.transpose()
    except FileNotFoundError:
        print('Warning: no cached codebooks found')

    # the orphan data
    model_results.orphan_labels = pd.DataFrame(
        data=orphan_label_array[:, :13], columns=test_label_columns
    )

    model_results.orphan_data = orphan_data
    model_results.orphan_vq2 = orphan_vq2
    model_results.orphan_vq1_ind = orphan_vq1_ind
    model_results.orphan_vq2_ind = orphan_vq2_ind
    return model_results
