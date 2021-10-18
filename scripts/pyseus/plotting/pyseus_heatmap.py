import sys
sys.path.append('../')
import matplotlib.pyplot as plt
import matplotlib
from numbers import Number
import numpy as np
import pandas as pd
from pyseus import basic_processing as pys
import plotly.offline
from plotly import graph_objs as go
import seaborn as sns
import plotly.figure_factory as ff
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.cluster import KMeans
import time
import pdb


def subtract_prey_median(imputed_df, mad_mod=True, mad_factor=1):
    """As an option to visualize clustering so that each intensity
    is subtracted by the prey group median, this function
    alters the base dataframe with the transformation"""

    transformed = imputed_df.copy()


    # Remove info columns for now, will add back again after transformation
    transformed.drop(columns=['Info'], level='Baits', inplace=True)
    transformed = transformed.T

    # Get a list of the columns (baits or preys)
    cols = list(transformed)

    # go through each prey (now in columns) and subtract median
    for col in cols:
        transformed[col] = transformed[col] - transformed[col].median()
        if mad_mod:
            mad = transformed[col].mad() * mad_factor
            transformed[col] = transformed[col].apply(lambda x: x if x > mad else 0)

    transformed = transformed.T
    # if mad_mod:
    #     t_cols = list(transformed)
    #     for col in t_cols:
    #         mad = transformed[col].mad() * mad_factor
    #         transformed[col] = transformed[col].apply(lambda x: x if x > mad else 0)

    # transpose back to original shape and add the info columns again
    info_cols = list([col for col in list(imputed_df) if col[0] == 'Info'])
    for col in info_cols:
        transformed[col] = imputed_df[col]

    return transformed


def prey_kmeans(imputed_df, k=20, method='single', ordering=True, verbose=True):
    """Create a large k clustered groups, and sort them by average group intensity.
    Return a list of Protein IDs after the sort

    rtype: dendro_side plotly figurefactory
    rtype: dendro_leaves list"""

    if verbose:
        print("Generating prey hierarchies and dendrogram...")
        start_time = time.time()

    # Create a median_df, taking median of all replicates
    median_df = pys.median_replicates(imputed_df, save_info=True, col_str='')
    median_df.drop(columns=['Protein names', 'Gene names',
    'Majority protein IDs'], inplace=True)

    # Protein IDs will be the reference to retrieve the correct order of preys
    median_df.set_index('Protein IDs', inplace=True)

    # Conduct K means clustering
    kmeans_model = KMeans(n_clusters=k).fit(median_df)
    kmeans_clusters = kmeans_model.predict(median_df)

    median_df['cluster'] = kmeans_clusters

    # Sort clusters by cluster average intensity
    grouped_df = median_df.groupby(['cluster'])
    cluster_intensities = grouped_df.mean()

    # Create a hierarchy of the clusters
    cluster_linkage = linkage(cluster_intensities, method=method,
        optimal_ordering=ordering)

    # List of clusters to be plotted sequentially
    cluster_leaves = leaves_list(cluster_linkage)

    # list of preys to be populated from cluster sequence
    leaves = []

    # sort thrugh clusters and populate with hierarchy of individual leaves
    for cluster in cluster_leaves:
        cluster_df = median_df[median_df['cluster'] == cluster]
        cluster_df.drop(columns=['cluster'], inplace=True)

        if cluster_df.shape[0] > 1:
            # Use plotly function to generate a linkage
            prey_linkage = linkage(cluster_df, method=method, optimal_ordering=ordering)

            # Retrieve the order of preys in the new linkage
            prey_leaves = leaves_list(prey_linkage)
            prey_leaves = [list(cluster_df.index)[x] for x in prey_leaves]

        else:
            prey_leaves = list(cluster_df.index)

        # add to the master list of leaves
        leaves = leaves + prey_leaves

    if verbose:
        end_time = np.round(time.time() - start_time, 2)
        print("Finished generating linkage in " + str(end_time) + " seconds.")

    return leaves


def bait_leaves(imputed_df, method='average', distance='euclidean', verbose=True):
    """Calculate the prey linkage and return the list of
    prey plotting sequence to use for heatmap. Use prey_kmeans for better performance
    rtype: prey_leaves list"""

    if verbose:
        print("Generating bait linkage...")
        start_time = time.time()
    # Create a median_df, taking median of all replicates
    median_df = pys.median_replicates(imputed_df, save_info=True, col_str='')
    median_df.drop(columns=['Protein names', 'Gene names',
    'Protein IDs', 'Majority protein IDs'],
        inplace=True)

    # Transpose to get linkages of baits
    median_df = median_df.T

    bait_linkage = linkage(median_df, method=method, optimal_ordering=True)

    # Retreieve the order of baits in the new linkage
    bait_leaves = leaves_list(bait_linkage)
    bait_leaves = [list(median_df.index)[x] for x in bait_leaves]

    if verbose:
        end_time = np.round(time.time() - start_time, 2)
        print("Finished generating linkage in " + str(end_time) + " seconds.")

    return bait_leaves


def prey_leaves(imputed_df, method='average', distance='euclidean', verbose=True):
    """Calculate the prey linkage and return the list of
    prey plotting sequence to use for heatmap. Use prey_kmeans for better performance.

    rtype: prey_leaves list"""
    if verbose:
        print("Generating prey linkage...")
        start_time = time.time()

    # Create a median_df, taking median of all replicates
    median_df = pys.median_replicates(imputed_df, save_info=True, col_str='')
    median_df.drop(columns=['Fasta headers', 'Protein names', 'Gene names'],
        inplace=True)

    # Protein IDs will be the reference to retrieve the correct order of preys
    median_df.set_index('Protein IDs', inplace=True)


    prey_linkage = linkage(median_df, method=method)

    # Retrieve the order of preys in the new linkage
    prey_leaves = leaves_list(prey_linkage)
    prey_leaves = [list(median_df.index)[x] for x in prey_leaves]


    if verbose:
        end_time = np.round(time.time() - start_time, 2)
        print("Finished generating linkage in " + str(end_time) + " seconds.")

    return prey_leaves


def dendro_heatmap(imputed_df, prey_leaves, hexmap, zmin, zmax, bait_leaves=None,
        bait_clust=False, verbose=True):
    """ From the dendro_leaves data, generate a properly oriented
    heatmap

    rtype fig pyplot Fig"""

    if verbose:
        print("Generating Heatmap...")
        start_time = time.time()

    plot_df = imputed_df.copy()

    # Set index to Protein IDs to match the dendro leaves
    plot_df.set_index(('Info', 'Protein IDs'), inplace=True)

    # Correctly order the plot df according to dendro leaves
    plot_df = plot_df.T[prey_leaves].T

    # Reset index to Gene Names
    plot_df.set_index(('Info', 'Gene names'), inplace=True)

    # Reorder columns based on bait_leaves
    if bait_clust:
        plot_df = plot_df[bait_leaves]

    # Informational columns are unnecessary now, drop them
    plot_df.drop(columns=['Info'], level='Baits', inplace=True)
    plot_df = plot_df.droplevel('Baits', axis=1)
    plot_df.index.rename('Gene names', inplace=True)



    # Generate the heatmap
    heatmap = [
        go.Heatmap(x=list(plot_df), y=list(plot_df.index), z=plot_df.values.tolist(),
        colorscale=hexmap, zmin=zmin, zmax=zmax)]



    if verbose:
        end_time = np.round(time.time() - start_time, 2)
        print("Finished heatmap in " + str(end_time) + " seconds.")

    return heatmap


def df_min_max(df):
    """Quickly output min and max values of the df"""

    # flatten the df to a list of all values
    all_vals = df.values.flatten().tolist()
    all_vals = list(filter(lambda x: isinstance(x, Number), all_vals))

    return min(all_vals), max(all_vals)


def color_map(df, zmin, zmax):
    """generate a color map, zmin, and zmax that the heatmap function will use
    Will add customization features in the future"""

    dfmin, dfmax = df_min_max(df)
    if zmin is None:
        zmin = dfmin
    if zmax is None:
        zmax = dfmax

    # Use built in seaborn function to blend palette
    cmap = sns.blend_palette(('black', 'blue', 'green', 'yellow',
        'orange', 'red'), n_colors=8, as_cmap=False)
    hexmap = []
    # change to hex values that Plotly can read
    for color in cmap:
        hexmap.append(matplotlib.colors.rgb2hex(color))

    # a range list from zmin to zmax
    a = list(range(int(zmin), int(zmax)))
    y = [0]*len(a)
    # plot colorscale
    fig = go.Figure(go.Heatmap(x=a, y=y, z=a, zmin=zmin, zmax=zmax,
        colorscale=hexmap, showscale=False),
        layout=go.Layout({'width': 1000, 'height': 200}, yaxis={'showticklabels': False}))

    return fig, zmin, zmax, hexmap



    # Keep only numerics


#     return heatmap_fig

# def heatmap_fig(df,hexmap,cols,width=600,height=800):
#     """generate a plotly heatmap given the df and selected columns"""
#     fig = go.Figure(data=go.Heatmap(
#                 z= df[cols].values.tolist(),
#                 x= cols,
#                 y= list(df.index),
#                 colorscale = hexmap,
#                 reversescale=False,
#                 hoverongaps = False,
#                 zmin=5,
#                 zmax= 14),
#                layout = go.Layout(width = width,height = height))
#     fig.show()
