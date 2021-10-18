import pandas as pd
import networkx as nx
import numpy as np
from itertools import repeat
from multiprocessing import Pool

import markov_clustering as mc


def retrieve_cluster_df(network, cluster, target_col='target', prey_col='prey'):
    """
    From a  list of cluster members from clusterone results, retrieve
    corresponding cluster from the network
    """

    network = network.copy()
    network = network[network[target_col] != network[prey_col]]
    return network[(network[target_col].isin(cluster)) & (network[prey_col].isin(cluster))]


def mcl_haircut(first_mcl, all_hits, target_col, prey_col, edge='', edge_thresh=0, clean=True):
    """
    Haircut clusterone members (removing single edge interactors)
    """
    all_hits = all_hits.copy()
    first_mcl = first_mcl.copy()
    members = first_mcl['gene_name']

    multi_args = zip(
        members, 
        repeat(all_hits), 
        repeat(target_col), 
        repeat(prey_col), 
        repeat(edge), 
        repeat(edge_thresh)
    )

    # multi processing for recursive haircut
    p = Pool()
    haircuts = p.starmap(recursive_haircut, multi_args)
    p.close()
    p.join()

    first_mcl['haircut_members'] = haircuts

    if clean:
        haircut_final = first_mcl[first_mcl['haircut_members'].apply(len) >= 2].reset_index(drop=True)
        haircut_final.drop(columns=['gene_name', 'mcl_cluster'], inplace=True)
        haircut_final.rename(columns={'haircut_members':'gene_names'}, inplace=True)
        haircut_final.reset_index(inplace=True)
        haircut_final.rename(columns={'index': 'super_cluster'}, inplace=True)

        return haircut_final
    
    else:
        return first_mcl


def recursive_haircut(cluster, all_hits, target_col, prey_col, edge='', edge_thresh=0):
    """
    Continue haircuts until there are no more single edges left
    """

    orig_len = len(cluster)
    new_len = 0
    clipped = cluster

    # while the haircut function clips more members from the cluster, continue 
    # applying the haircut function
    while orig_len - new_len > 0:
        orig_len = len(clipped)
        clipped = haircut(clipped, all_hits, target_col, prey_col, edge, edge_thresh)
        new_len = len(clipped)

    return clipped


def haircut(cluster, all_hits, target_col, prey_col, edge, edge_thresh):
    """
    a single iteration of a haircut to remove single edge interactors
    """

    # Find all edges that belong to a given cluster
    cluster_group = all_hits[
        (all_hits[target_col].isin(cluster)) & (all_hits[prey_col].isin(cluster))
    ]

    haircut_members = []
    for member in cluster:
        membership = cluster_group[
            (cluster_group[target_col] == member) | (cluster_group[prey_col] == member)
        ]

        # validation that a cluster member is involved in at least two interactions
        members = set(membership[target_col].to_list() + membership[prey_col].to_list())
        
        # if there are more than two edges, retain the cluster member
        if len(members) > 2:
            haircut_members.append(member)

        # optionally, if the edge weight of a single edge is higher than 
        # the given threshold, retain the cluster member
        elif edge:
            if membership[edge].max() > edge_thresh:
                haircut_members.append(member)

    return haircut_members


def clean_up_cytoscape_mcl(mcl_stoi, grouped=False):
    """
    cytoscape MCL output has quirky column names- rename to mcl_cluster
    and gene_name columns. Also remove clusters that have less than 2 members.
    Group option is available to group protein members in a list by cluster 
    """

    mcl_stoi = mcl_stoi.copy()
    
    # drop unnecessary columns and rename column names
    mcl_stoi.drop(columns=['selected', 'shared name'], inplace=True)
    mcl_stoi.reset_index(drop=True, inplace=True)
    mcl_stoi.rename(columns={'__mclCluster': 'mcl_cluster', 'name': 'gene_name'}, inplace=True)
    mcl_stoi.sort_values(['mcl_cluster', 'gene_name'], inplace=True)
    mcl_stoi = mcl_stoi[mcl_stoi['mcl_cluster'].apply(lambda x: np.isfinite(x))]

    # drop clusters that have less than 2 members
    count = mcl_stoi.groupby('mcl_cluster').count()
    clusters = count[count['gene_name'] > 2].index.to_list()
    mcl_stoi = mcl_stoi[mcl_stoi['mcl_cluster'].isin(clusters)]
    mcl_stoi.reset_index(drop=True, inplace=True)
    
    # grouping option 
    if grouped:
        mcl_stoi = pd.DataFrame(
            mcl_stoi.groupby('mcl_cluster')['gene_name'].apply(list)).reset_index()
    
    return mcl_stoi


def second_mcl(
    first_mcl, 
    network, 
    target_col, 
    prey_col, 
    first_thresh, 
    mcl_thresh, 
    mcl_inflation, 
    edge='', 
    clean=True
):
    """
    Performs secondary clustering from the first MCL cluster results. 
    Require network df that contains PPI edges and the cleaned first_mcl results
    """
    network = network.copy()
    first_mcl = first_mcl.copy()
    clusters = first_mcl['gene_names']

    # fist clustering id
    c_clusters = []

    # Boolean for whether the cluster went thru second clustering
    mcl = []

    # New MCL cluster
    m_clusters = []

    for idx, cluster in enumerate(clusters):

        # go thru the MCL second cluster
        cluster_network = retrieve_cluster_df(network, cluster, target_col, prey_col)
        
        # if a cluster has less than 2 interactions, do not cluster again
        if cluster_network.shape[0] < 2:
            c_clusters.append(idx)
            mcl.append(False)
            m_clusters.append([])
            continue
        else: 
            if not edge:
                edge = None

            # NetworkX transformation of pandas interactions to sparse matrix
            c_graph = nx.convert_matrix.from_pandas_edgelist(
                cluster_network, target_col, prey_col, edge_attr=edge
            )
            nodes = list(c_graph.nodes)
            c_mat = nx.to_scipy_sparse_matrix(c_graph)

            # Run MCL with a given inflation parameter
            result = mc.run_mcl(c_mat, inflation=mcl_inflation)
            mcl_clusters = mc.get_clusters(result, keep_overlap=False)

        # append new second_clustering members
        for mcl_cluster in mcl_clusters:
            if len(mcl_cluster) >= mcl_thresh:
                mcl_nodes = [nodes[x] for x in mcl_cluster]
                c_clusters.append(idx)
                mcl.append(True)
                m_clusters.append(mcl_nodes)

    # organize second clustering results into a dataframe
    mcl_df = pd.DataFrame()
    mcl_df['super_cluster'] = c_clusters
    mcl_df['second_clustering'] = mcl
    mcl_df['gene_names'] = m_clusters

    # clean up the table to exclude empty core clusters and explode gene names
    if clean:
        mcl_df = mcl_df[mcl_df['gene_names'].apply(len) > 0].reset_index(drop=True)
        core = mcl_df.explode('gene_names').reset_index().rename(columns={'index': 'core_cluster'})
        core.drop(columns=['second_clustering'], inplace=True)

        return core
    else:
        return mcl_df


def split_mcl(
    summary, 
    network, 
    target_col, 
    prey_col, 
    mcl_thresh, 
    mcl_inflation,
    edge='', 
    raw_return=False
):
    """
    If MCL cluster has two core clusters, try to split it into two
    """
    network = network.copy()
    summary = summary.copy()

    clusters = summary.groupby('super_cluster')['gene_name'].apply(list).to_list()
    num_cores = pd.DataFrame(summary.groupby('super_cluster')['core_cluster'].nunique()).reset_index()
    change_list = num_cores[num_cores['core_cluster'] >= 2]['super_cluster'].to_list()

    # first clustering id
    c_clusters = []

    # New MCL cluster
    m_clusters = []
    for idx, cluster in enumerate(clusters):

        # original clusterone cluster number
        idx += 1
        if idx not in change_list:
            m_clusters.append(cluster)
            c_clusters.append(idx)
            continue

        # go thru the MCL second cluster
        else:
            cluster_network = retrieve_cluster_df(network, cluster, target_col, prey_col)
            if not edge:
                edge = None
    
            # NetworkX transformation of pandas interactions to sparse matrix
            c_graph = nx.convert_matrix.from_pandas_edgelist(
                cluster_network, target_col, prey_col, edge_attr=edge
            )
            nodes = list(c_graph.nodes)
            c_mat = nx.to_scipy_sparse_matrix(c_graph)

            # Run MCL with a given inflation parameter
            result = mc.run_mcl(c_mat, inflation=mcl_inflation)
            mcl_clusters = mc.get_clusters(result, keep_overlap=False)

            for mcl_cluster in mcl_clusters:
                if len(mcl_cluster) >= mcl_thresh:
                    mcl_nodes = [nodes[x] for x in mcl_cluster]
                    c_clusters.append(idx)
                    m_clusters.append(mcl_nodes)

    mcl_df = pd.DataFrame()
    # mcl_df['super_cluster'] = c_clusters
    mcl_df['gene_name'] = m_clusters
    mcl_df.reset_index(inplace=True)
    mcl_df.rename(columns={'index':'super_cluster'}, inplace=True)

    if raw_return:
        return mcl_df

    new = mcl_df.explode('gene_name').reset_index().rename(
        columns={'index': 'renewed_super_cluster'}
    )
    new['renewed_super_cluster'] = new['renewed_super_cluster'].apply(lambda x: x+1)
    masters = summary.merge(new, on=['super_cluster', 'gene_name'], how='left')
    masters = masters.sort_values(
        ['renewed_super_cluster', 'core_cluster', 'gene_name']
    ).reset_index(drop=True)
    
    return masters
