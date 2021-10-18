import os
import umap
import dask
import datetime
import skimage
import numpy as np
import scipy as sp
import pandas as pd
import anndata as ad

import sklearn
import sklearn.manifold
import sklearn.decomposition

import dask.diagnostics
import plotly.graph_objects as go
import plotly.express as px

from scipy.cluster.hierarchy import dendrogram
from matplotlib import pyplot as plt

from . import analysis_utils

SCATTERPLOT_COLORS = [
    '#f44336', '#e91e63', '#9c27b0', '#673ab7', '#3f51b5',
    '#2196f3', '#03a9f4', '#00bcd4', '#009688', '#4caf50',
    '#8bc34a', '#cddc39', '#ffeb3b', '#ffc107', '#ff9800',
    '#ff5722', '#795548', '#009e9e', '#607d8b', '#d500f9',
    '#212121', '#ff9e80', '#ff6d00', '#ffff00', '#76ff03',
    '#00e676', '#64ffda', '#18ffff',
]


def timestamp():
    return datetime.datetime.now().strftime('%Y-%m-%d')


class ModelResults:

    def __init__(
        self, 
        nickname,
        num_features,
        test_data,
        test_labels,
        test_vq1,
        test_vq2,
        test_vq1_ind,
        test_vq2_ind,
        exports_dirpath=None
    ):

        # dataset nickname (used for generate figure filenames)
        self.nickname = nickname

        # the number of features in the codebook
        self.num_features = num_features

        # where to save final/processed results
        self.exports_dirpath = exports_dirpath
        if self.exports_dirpath:
            os.makedirs(self.exports_dirpath, exist_ok=True)

        self.test_data = test_data
        self.test_labels = test_labels
        self.test_vq1 = test_vq1
        self.test_vq2 = test_vq2
        self.test_vq1_ind = test_vq1_ind
        self.test_vq2_ind = test_vq2_ind

        self.target_histograms = {}
        self.target_vectors = {}

        # categorical for scatterplots
        self.scatterplot_colors = SCATTERPLOT_COLORS
        self.scatterplot_colors = np.array(self.scatterplot_colors)[::2]

    
    @staticmethod
    def construct_labels_dataframe(labels_array):
        '''
        labels_array: numpy array of data (patch) labels
        '''
        labels = pd.DataFrame(
            data=labels_array, 
            columns=[
                'cell_line_id', 'target_name', 'unknown_id', 'filepath',
                'label_0', 'label_1', 'label_2', 'label_3', 'label_4', 'label_5',
            ]
        )
        return labels


    def construct_target_labels(self):
        '''
        Generate a dataframe of target labels from the dataframe of test (patch) labels
        '''
        target_labels = self.test_labels.groupby('cell_line_id').first().reset_index()

        # drop patch-specific columns
        target_labels.drop(
            labels=['pml_id', 'fov_id', 'umap_0', 'umap_1', 'filepath'], 
            axis=1, 
            inplace=True,
            errors='ignore'
        )
        return target_labels


    def concatenate_orphans(self):
        '''
        Concatenate the orphan test labels and VQ2 vectors
        NOTE: this is specific to december-2020 results
        '''

        assert self.orphan_vq2.shape[0] == self.orphan_vq2_ind.shape[0]
        num_orphans = self.orphan_vq2.shape[0]

        self.test_vq2 = np.concatenate((self.test_vq2, self.orphan_vq2), axis=0)
        self.test_vq1_ind = np.concatenate((self.test_vq1_ind, self.orphan_vq1_ind), axis=0)
        self.test_vq2_ind = np.concatenate((self.test_vq2_ind, self.orphan_vq2_ind), axis=0)
        
        self.test_labels = pd.concat(
            (self.test_labels, self.orphan_labels.iloc[:num_orphans]),
            axis=0
        )

        if self.test_data is not None:
            # the test data has three channels, but the orphan data has two
            self.test_data = np.concatenate((self.test_data[:, :, :, :2], self.orphan_data), axis=0)


    def figure_filepath(self, tag, vq=None):
        '''
        '''
        nickname = '%s-%s' % (self.nickname, vq) if vq is not None else self.nickname
        return os.path.join('..', 'figures', '%s_%s' % (nickname, tag))


    def _construct_codebook(self, vq):
        '''
        '''
        test_vq = self.test_vq1 if vq == 'vq1' else self.test_vq2
        test_vq_ind = self.test_vq1_ind if vq == 'vq1' else self.test_vq2_ind

        num_features = self.num_features
        num_feature_dims = test_vq.shape[-1]

        codebook = np.zeros((num_features, num_feature_dims))
        tmp_inds = test_vq_ind[:10000, :, :]
        tmp_vecs = test_vq[:10000, :, :, :]

        feature_inds = np.unique(tmp_inds)
        for ind in feature_inds:
            pos = np.argwhere(tmp_inds == ind)[0]
            codebook[ind, :] = tmp_vecs[pos[0], pos[1], pos[2], :]
            
        return codebook


    def construct_codebooks(self):
        '''
        '''
        self.codebook_vq1 = self._construct_codebook(vq='vq1')
        self.codebook_vq2 = self._construct_codebook(vq='vq2')


    def get_codebook_nonzero(self, vq):
        codebook = self.codebook_vq1 if vq == 'vq1' else self.codebook_vq2
        return codebook[np.abs(codebook).sum(axis=1) != 0, :]


    def plot_codebook(self, vq, figsize=(6, 6), save=False):
        '''
        '''
        codebook = self.codebook_vq1 if vq == 'vq1' else self.codebook_vq2
        codebook_nonzero = self.get_codebook_nonzero(vq)
        maxx = np.abs(codebook).max()

        plt.figure(figsize=figsize)
        plt.imshow(codebook, cmap='RdBu', vmin=-maxx, vmax=maxx)
        if save:
            plt.savefig(self.figure_filepath(vq=vq, tag='codebook-all.pdf'))

        plt.figure(figsize=figsize)
        plt.imshow(codebook_nonzero, cmap='RdBu', vmin=-maxx, vmax=maxx)
        if save:
            plt.savefig(self.figure_filepath(vq=vq, tag='codebook-nonzero.pdf'))


    def simplify_codebook(self, vq, num_clusters):
        '''
        num_clusters : the number of clusters into which to cluster/group the features
        '''

        codebook = self.codebook_vq1 if vq == 'vq1' else self.codebook_vq2
        num_features = codebook.shape[0]
        model = sklearn.cluster.AgglomerativeClustering(
            distance_threshold=None, 
            n_clusters=num_clusters
        )
        model = model.fit(codebook)
        return model.labels_
    

    def plot_codebook_dendrogram(self, vq):
        '''
        '''

        codebook_nonzero = self.get_codebook_nonzero(vq)
        model = sklearn.cluster.AgglomerativeClustering(distance_threshold=0, n_clusters=None)
        model = model.fit(codebook_nonzero)

        # the block below was copied from sklearn docs
        # https://scikit-learn.org/stable/auto_examples/cluster/plot_agglomerative_dendrogram.html    
    
        # create the counts of samples under each node
        counts = np.zeros(model.children_.shape[0])
        n_samples = len(model.labels_)
        for i, merge in enumerate(model.children_):
            current_count = 0
            for child_idx in merge:
                if child_idx < n_samples:
                    current_count += 1  # leaf node
                else:
                    current_count += counts[child_idx - n_samples]
            counts[i] = current_count
        linkage_matrix = np.column_stack([model.children_, model.distances_, counts])

        plt.figure(figsize=(6, 15))
        plt.rcParams.update({'font.size': 12})
        d = dendrogram(
            linkage_matrix, 
            orientation='left', 
            leaf_font_size=10, 
            truncate_mode=None
        )

        # show the code book with rows in the order of the dendrogram 
        inds = np.array(d['ivl']).astype(int)
        plt.figure(figsize=(10, 10))
        plt.imshow(codebook_nonzero[inds, :])
        return d


    def tile_patches(self, target_name, vq, offset=0, background_inds=None):
        
        patch_size = 100
        test_vq_ind = self.test_vq1_ind if vq == 'vq1' else self.test_vq2_ind
        bottleneck_size = test_vq_ind.shape[-1]
        
        n_rows, n_cols = 3, 5
        inds_tile = np.zeros((n_rows*bottleneck_size, n_cols*bottleneck_size), dtype='uint8')
        patch_tile = np.zeros((n_rows*patch_size, n_cols*patch_size), dtype='uint8')

        mask = self.test_labels.target_name == target_name
        masked_data = self.test_data[mask, :]
        masked_inds = test_vq_ind[mask, :]

        count = 0
        for row in range(n_rows):
            for col in range(n_cols):
                inds_tile[
                    row*bottleneck_size:(row + 1)*bottleneck_size, 
                    col*bottleneck_size:(col + 1)*bottleneck_size
                ] = masked_inds[count, :, :]
                
                patch_tile[
                    row*patch_size:(row + 1)*patch_size, 
                    col*patch_size:(col + 1)*patch_size
                ] = analysis_utils.autoscale(masked_data[count, :, :, 0], percentile=1)
                
                count += 1
        
        # upsample the indices tile
        inds_tile = skimage.transform.resize(
            inds_tile, patch_tile.shape, order=0, preserve_range=True
        )

        background_tile = (inds_tile * 0).astype(bool)
        for ind in background_inds:
            background_tile = (inds_tile == ind) | background_tile

        return patch_tile, inds_tile, background_tile


    def construct_target_histograms(self, vq):
        '''
        '''
        test_vq_ind = self.test_vq1_ind if vq == 'vq1' else self.test_vq2_ind

        # with these edges, feature counts are in counts[feature_ind],
        # assuming that the feature indices are in range(0, num_features)
        bin_edges = np.arange(0, self.num_features + 1)

        target_labels = self.construct_target_labels()
        counts = np.zeros((target_labels.shape[0], len(bin_edges) - 1))
        for ind, row in target_labels.iterrows():
            mask = self.test_labels.cell_line_id == row.cell_line_id
            _counts, _edges = np.histogram(
                test_vq_ind[mask, :].flatten(), bins=bin_edges, density=True
            )
            counts[ind, :] = _counts

        self.target_histograms[vq] = {'counts': counts, 'target_labels': target_labels}


    def construct_target_vectors(self, vq, n_components=50, using='mean'):
        '''
        Construct a per-target vector directly from the VQ vectors themselves
        n_components : the number of PCs to average for the per-target vector;
            if None, no PCA is performed (and the full 9216-dimensional vector is averaged)
            empirically, 50-100 for a VQ2 bottleneck of 12x12x64 seems to be enough
        '''
        X = self.test_vq1 if vq == 'vq1' else self.test_vq2
        
        # flatten over the bottleneck dimensions (usually something like 12x12x64)
        X = X.reshape(X.shape[0], -1)

        # take the first n PCs to lower the dimensionality
        # at least to something less than the number of targets
        if n_components is not None:
            model = sklearn.decomposition.PCA(n_components=n_components)
            X_reduced = model.fit_transform(X)
        else:
            X_reduced = X

        target_labels = self.construct_target_labels()
        averaged_vectors = np.zeros((target_labels.shape[0], X_reduced.shape[1]))
        for ind, row in target_labels.iterrows():
            mask = self.test_labels.cell_line_id == row.cell_line_id
            if using == 'mean':
                averaged_vectors[ind, :] = np.mean(X_reduced[mask, :], axis=0)
            elif using == 'median':
                averaged_vectors[ind, :] = np.median(X_reduced[mask, :], axis=0)

        self.target_vectors[vq] = {'X': averaged_vectors, 'target_labels': target_labels}


    @staticmethod
    def plot_explained_variance(X, n_dims=None):
        model = sklearn.decomposition.PCA(n_components=(n_dims or min(X.shape)))
        coords = model.fit_transform(X) 
        plt.plot(model.explained_variance_ratio_)
        return model
    

    def embed_target_histograms(
        self, 
        vq, 
        num_codebook_clusters=None, 
        num_pca_dims=None, 
        **umap_kwargs
    ):
        '''
        num_codebook_clusters : the number of clusters into which to condense the codebook features
            and to aggregate the histograms
            Note: no clustering/condensation if None
        num_pca_dims : the number of PCA components to retain 
            from the (possibly codebook-clustering-condensed) raw histograms
            Note: no PCA if None
        '''

        tags = ['coords']
        
        counts = self.target_histograms[vq]['counts'].copy()

        if num_codebook_clusters is not None:
            labels = self.simplify_codebook(vq=vq, num_clusters=num_codebook_clusters)
            counts_reduced = np.zeros((counts.shape[0], max(labels) + 1))
            for label in labels:
                counts_reduced[:, label] = counts[:, labels == label].sum(axis=1)
            counts = counts_reduced
            tags.append('cc%d' % num_codebook_clusters)

        if num_pca_dims is not None:
            model = sklearn.decomposition.PCA(n_components=num_pca_dims)
            coords = model.fit_transform(counts)
            tags.append('pc%d' % num_pca_dims)
        else:
            coords = counts

        name = '_'.join(tags)
        self.target_histograms[vq][name] = coords

        default_umap_kwargs = dict(n_neighbors=10, min_dist=0.1, spread=1)
        default_umap_kwargs.update(**umap_kwargs)
        reducer = umap.UMAP(**default_umap_kwargs)
        reducer.fit(coords)
        if umap_kwargs:
            tags.append('umap-custom')
        else:
            tags.append('umap')

        name = '_'.join(tags)
        self.target_histograms[vq][name] = reducer.embedding_.copy()


    def get_positions(self, vq, cleanup_method, show_pca=False):
        '''
        Retrieve the specified 'kind' of embedded target histograms
        (by PCA or by umap)

        cleanup_method : whether the input histograms where cleaned up by 
            PCA or codebook clustering
        '''
        # return the PCA coordinates themselves
        if show_pca:
            positions = self.target_histograms[vq]['pca_coords'].copy()*15

        # return the umap coordinates from the specified kind of cleaned-up histograms
        else:
            positions = self.target_histograms[vq]['%s_umap_reducer' % cleanup_method].embedding_.copy()

        target_labels = self.target_histograms[vq]['target_labels'].copy()

        return positions, target_labels


    def plot_embedding(
        self, vq, name=None, user_provided_positions=None, ax=None, legend=True, **scatter_kwargs
    ):
        '''
        Create a 2D scatterplot of targets using either PCA or umap coordinates
        '''
        if ax is None:
            plt.figure(figsize=(12, 12))
            plt.rcParams.update({'font.size': 12})
            ax = plt.gca()

        target_labels = self.target_histograms[vq]['target_labels'].copy()
        target_names = target_labels.target_name

        if user_provided_positions is not None:
            positions = user_provided_positions.copy()
        else:
            positions = self.target_histograms[vq][name].copy()

        # find the most frequent label_0 values
        label_counts = target_labels.label_0.value_counts()
        most_common_labels = label_counts[:len(self.scatterplot_colors)].index

        for ind, label in enumerate(most_common_labels):
            target_names_with_label = target_labels.loc[target_labels.label_0 == label].target_name.values
            positions_mask = target_names.isin(target_names_with_label)
            ax.scatter(
                positions[positions_mask, 0], 
                positions[positions_mask, 1], 
                color=self.scatterplot_colors[ind], 
                label=label,
                **scatter_kwargs
            )
        if legend:
            ax.legend()

        ax.set_title(name if name else 'user provided positions')


    def plot_embedding_with_patches(
        self, vq, cleanup_method='clustering', show_pca=False, user_provided_positions=None, min_num_patches=0
    ):
        '''
        '''
        positions, target_labels = self.get_positions(vq, cleanup_method, show_pca)
        target_names = target_labels.target_name

        if user_provided_positions is not None:
            positions = user_provided_positions.copy()

        reference_patches = []
        for target_name in target_names:
            ind = self.test_labels.loc[self.test_labels.target_name == target_name].index[0]
            patch = analysis_utils.autoscale(self.test_data[ind, :, :, 0], percentile=1)
            reference_patches.append(patch)

        fig = plt.figure(figsize=(12, 12))
        for pos, patch, target_name in zip(positions, reference_patches, target_names):
            ax = fig.add_axes([pos[1], pos[0], 1/8, 1/8], frameon=False)
            ax.imshow(patch[::2, ::2], cmap='gray')
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(target_name, fontsize=18)


    def plot_embedding_with_plotly(
        self, vq, cleanup_method='clustering', show_pca=False, user_provided_positions=None
    ):
        '''
        '''
        positions, target_labels = self.get_positions(vq, cleanup_method, show_pca)
        target_names = target_labels.target_name

        if user_provided_positions is not None:
            positions = user_provided_positions.copy()

        labels = self.test_labels.copy()

        # array of color indices for each target
        colors = np.zeros((len(target_names),))

        # the label_0 value for each target
        target_first_labels = np.array(['unknown']*len(target_names)).astype('<U20')

        for ind, label in enumerate(labels.label_0.unique()):
            target_names_with_label = labels.loc[labels.label_0 == label].target_name.values
            mask = target_names.isin(target_names_with_label)
            
            colors[mask] = ind + 1
            target_first_labels[mask] = label

        # concatenate all of the labels for each target
        labels_grouped = labels.groupby('target_name').first()
        label_columns = ['label_0', 'label_1', 'label_2', 'label_3', 'label_4', 'label_5']
        target_all_labels = []
        for target_name in target_names:
            all_labels = labels_grouped.loc[target_name][label_columns].values
            target_all_labels.append(', '.join(all_labels[~pd.isna(all_labels)]))

        df = pd.DataFrame(
            data={
                'x': positions[:, 0], 
                'y': positions[:, 1], 
                'name': target_names,
                'first_label': target_first_labels, 
                'all_labels': target_all_labels
            }
        )
        fig = px.scatter(
            df,
            x='x',
            y='y',
            color='first_label',
            color_discrete_sequence=px.colors.qualitative.D3,
            hover_name='name', 
            hover_data={'all_labels': True, 'x': False, 'y': False,},
        )
        axes_props = dict(
            showline=True, 
            linewidth=1, 
            linecolor='#999', 
            mirror=True, 
            zeroline=True, 
            zerolinecolor='#eee', 
            gridcolor='#eee'
        )
        fig.update_xaxes(**axes_props)
        fig.update_yaxes(**axes_props)
        fig.update_layout(**dict(
                height=800, 
                width=1000, 
                plot_bgcolor='white', 
            )
        )
        fig.show()
        return fig


    def export_image_patches(self):
        '''
        Downsample the patches (the test_data) to uint8 for the umap-viewer
        '''
        tasks = [dask.delayed(analysis_utils.autoscale)(image[:, :, 0]) for image in self.test_data]
        with dask.diagnostics.ProgressBar():
            downsampled_images = dask.compute(*tasks)
        downsampled_images = np.array(downsampled_images)

        filepath = os.path.join(self.exports_dirpath, '%s-uint8-patch-images.npy' % timestamp())
        np.save(filepath, downsampled_images)
        print('Patch images saved to %s' % filepath)
        return downsampled_images


    def export_patch_labels(self):
        '''
        Save the dataframe of patch labels and umap coords for the umap viewer
        '''
        labels = self.test_labels.copy()
        filepath = os.path.join(self.exports_dirpath, '%s-patch-labels-umap.csv' % timestamp())      
        labels.to_csv(filepath)
        print('Patch labels saved to %s' % filepath)


    def get_target_labels_filepath(self):
        '''
        The filepath used in export_target_labels
        '''
        return os.path.join(self.exports_dirpath, 'target-labels-umap.csv')


    def export_target_labels(self, vq):
        '''
        Save the dataframe of target labels and umap coords for the umap viewer

        TODO: this is currently done in notebooks to enable easily choosing 
        the kind of umap coords and clustering results to append to the labels
        (and to show in the umap viewer)
        '''
        labels = self.target_histograms[vq]['target_labels'].copy()
        filepath = self.get_target_labels_filepath()
        return

    
    def export_adata(
        self, vq, kind='vectors', using='mean', pub_ready_only=False, rerun=False, filepath=None
    ):
        '''
        vq : 'vq1' or 'vq2'
        kind : 'vectors' or 'histograms'
        using : 'mean' or 'median'
        rerun : whether to re-run the target vector construction methods
        filepath : optional filepath to which to save the adata object
        '''
        
        if kind == 'vectors':
            if not self.target_vectors.get(vq) or rerun:
                self.construct_target_vectors(vq, n_components=None, using=using)
            X = self.target_vectors[vq]['X'].copy()
            target_labels = self.target_vectors[vq]['target_labels'].copy()

        elif kind == 'histograms':
            if not self.target_histograms.get(vq) or rerun:
                self.construct_target_histograms(vq)
            X = self.target_histograms[vq]['counts'].copy()
            target_labels = self.target_histograms[vq]['target_labels'].copy()

        if pub_ready_only:
            mask = target_labels.pub_ready
            X = X[mask, :]
            target_labels = target_labels.loc[mask]
            
        adata = ad.AnnData(
            sp.sparse.csc_matrix(X), 
            obs=target_labels, 
            var=pd.DataFrame(index=np.arange(X.shape[1])), 
            dtype='float'
        )

        if filepath:
            adata.write_h5ad(filepath)
        return adata
