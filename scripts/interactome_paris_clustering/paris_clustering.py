import pandas as pd
import requests


def calculate_intercluster_edges(cluster_membership, interactions):
    """
    For hierarchical clustering of existing clusters, Paris algorithm
    requires edge weights between clusters. 
    For this, we simply calculate existing # of interactions amongst members
    of two clusters.
    """

    interactions = interactions.copy()
    
    cluster_membership = cluster_membership.copy()

    # All the unique clusters in the PPI network, designated by Markov clustering
    clusters = cluster_membership['community'].unique()

    # Returned dataframe will have three columns - origin cluster, target cluster,
    # and sum of interactions between the two clusters
    clust_ones = []
    clust_twos = []
    inter_edges = []

    # iterate through clusters
    for clust_one in clusters:
        # second iteration through clusters
        for clust_two in clusters:
            # exclude overlapping cluster-cluster combinations
            if clust_one <= clust_two:
                continue

            clust_ones.append(clust_one)
            clust_twos.append(clust_two)

            # Find all genes in membership
            clust_one_genes = cluster_membership[
                cluster_membership['community'] == clust_one
            ]['gene_names'].to_list()

            clust_two_genes = cluster_membership[
                cluster_membership['community'] == clust_two
            ]['gene_names'].to_list()

            # Find all interactions between two clusters and sum the # interactions
            inter_1 = interactions[
                (interactions['prot_1'].isin(clust_one_genes)) &
                (interactions['prot_2'].isin(clust_two_genes))
            ].drop_duplicates()

            inter_2 = interactions[
                (interactions['prot_1'].isin(clust_two_genes)) & 
                (interactions['prot_2'].isin(clust_one_genes))
            ].drop_duplicates()

            total = inter_1.shape[0] + inter_2.shape[0]

            inter_edges.append(total)

    cluster_edges = pd.DataFrame()
    cluster_edges['cluster_1'] = clust_ones
    cluster_edges['cluster_2'] = clust_twos
    cluster_edges['intersection'] = inter_edges

    cluster_edges = cluster_edges[cluster_edges['intersection'] > 0]
    return cluster_edges


def query_panther(target_names, all_target_names, biological):
    """
    from a list of genes out of database of genes, calculate GO enrichment and FDR
    """

    url = 'http://pantherdb.org/services/oai/pantherdb/enrich/overrep'
    # dataset ids from http://pantherdb.org/services/oai/pantherdb/supportedannotdatasets

    panther_datasets = {
        "molecular_function": "ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF",
        "biological_process": "ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP",
        "cellular_component": "ANNOT_TYPE_ID_PANTHER_GO_SLIM_CC"
    }
    if biological == 'bp':
        panther_dataset = panther_datasets['biological_process']
    elif biological == 'cc':
        panther_dataset = panther_datasets['cellular_component']
    else:
        panther_dataset = panther_datasets['molecular_function']

    panther_human_organism_id = 9606
    params = {
        'geneInputList': (','.join(target_names)),
        'refInputList': (','.join(all_target_names)),
        'organism': panther_human_organism_id,
        'refOrganism': panther_human_organism_id,
        'annotDataSet': panther_dataset,
        'enrichmentTestType': 'FISHER',
        'correction': 'FDR'
    }

    result = requests.post(url, params)
    d = result.json()

    # these are the hits (already sorted by p-value)
    df = pd.DataFrame(data=d['results']['result'])
    df['go_term_id'] = df.term.apply(lambda s: s.get('id'))
    df['go_term_label'] = df.term.apply(lambda s: s.get('label'))
    df.drop('term', axis=1, inplace=True)
    return df
