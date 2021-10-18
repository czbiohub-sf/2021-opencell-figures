'''
This script generates the patch labels dataframe and the uint8-downsampled array of patch images
that are required for the umap viewer app
'''

import os
import scanpy as sc
from . import model_results, clustering_workflows, ground_truth_labels


# for december results, model is 'full' or 'no-channel-splitting'
res = model_results.ModelResults.load_december_results(
    root_dirpath=None, 
    dataset='full', 
    model='no-decoder',
    rep=3
)

# concatenate the orphan test labels, test data, and VQ2 vectors 
res.concatenate_orphans()

# merge the target annotations (with hard-coded path to cache of the opencell lines/ endpoint)
res.test_labels = ground_truth_labels.merge_opencell_annotations(
    df=res.test_labels, 
    filepath='/ML_group/KC/2021-04-15-opencell-lines-endpoint.json', 
    only_one=True
)

print('Exporting patch labels dataframe')
res.export_patch_labels()

print('Exporting downsampled patch images')
downsampled_images = res.export_image_patches()
print('Downsampled patch images shape: %s' % (downsampled_images.shape,))
