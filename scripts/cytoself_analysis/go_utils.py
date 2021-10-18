import json
import requests
import pandas as pd
from collections import namedtuple
from scipy import stats


def query_enrichr(gene_names):
    
    enrichr_url = 'http://maayanlab.cloud/Enrichr/addList'
    description = ''
    payload = {
        'list': (None, '\n'.join(gene_names)),
        'description': (None, description)
    }

    response = requests.post(enrichr_url, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)
    user_list_id = data['userListId']

    enrichr_url = 'http://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'

    # gene_set_library = 'GO_Cellular_Component_2018'
    # gene_set_library = 'GO_Biological_Process_2018'
    gene_set_library = 'Jensen_COMPARTMENTS'
    
    response = requests.get(
        enrichr_url + query_string % (user_list_id, gene_set_library)
     )
    if not response.ok:
        raise Exception('Error fetching enrichment results')
    data = json.loads(response.text)

    enrichr_results_columns = [
        'rank', 
        'go_term_name', 
        'p_value', 
        'z_score', 
        'combined_score', 
        'overlapping_genes', 
        'adjusted_p_value', 
        'old_p_value', 
        'old_adjusted_p_value'
    ]
    enrichr_results = pd.DataFrame(
        data=data[gene_set_library], columns=enrichr_results_columns
    )
    
    # count the number of genes with the label (for enrichr results only)
    enrichr_results['num_overlapping_genes'] = enrichr_results.overlapping_genes.apply(
        lambda s: len(s)
    )
    enrichr_results = enrichr_results.sort_values(by='p_value')
    return enrichr_results


def query_panther(target_names, reference_target_names, dataset_kind):
    '''
    Query the Panther API to obtain enriched GO terms for a set of gene names,
    relative to a set of reference gene names

    dataset : one of 'cc', 'bp', or 'mf' (see below)
    '''

    # dataset ids from http://pantherdb.org/services/oai/pantherdb/supportedannotdatasets
    panther_dataset_ids = {

        # molecular function
        'mf': 'GO:0003674',

        # biological process
        'bp': 'GO:0008150',

        # cellular component
        'cc': 'GO:0005575',
    }

    panther_human_organism_id = 9606
    url = 'http://pantherdb.org/services/oai/pantherdb/enrich/overrep'

    params = {
        'geneInputList': (','.join(target_names)),
        'refInputList': (','.join(reference_target_names)),
        'organism': panther_human_organism_id,
        'refOrganism': panther_human_organism_id,
        'annotDataSet': panther_dataset_ids[dataset_kind],
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


def calc_enrichment_pvals(df):
    '''
    Calculate enrichment of labels in clusters

    df: a dataframe of targets/genes with an ensg_id, cluster_id, and label_id column
        ensg_id : the ENSG ID of the target/gene
        cluster_id : arbitrary categorical variable that identifies the clusters
            in which to calculate the enrichment of the labels
            (we assume that each gene is in only one cluster)
        label_id : the labels (each gene may have one or more labels)
    '''
    label_counts = df.label_id.value_counts()

    # total number of targets with at least one label
    num_targets = len(df.ensg_id.unique())

    Result = namedtuple(
        'result', ['cluster_id', 'label_id', 'cluster_size', 'label_enrichment', 'corrected_pval']
    )

    rows = []
    for cluster_id in df.cluster_id.unique():
        dff = df.loc[df.cluster_id == cluster_id].copy()

        # labels in the cluster
        num_labels_in_cluster = len(dff.label_id.unique())

        # targets in the cluster
        num_targets_in_cluster = len(dff.ensg_id.unique())

        for label_id in dff.label_id.unique():
            
            num_targets_in_cluster_with_label = dff.label_id.value_counts()[label_id]

            # overall number of targets with the label
            num_targets_with_label = label_counts[label_id]

            table_a = num_targets_in_cluster_with_label
            table_b = num_targets_in_cluster - num_targets_in_cluster_with_label
            table_c = num_targets_with_label - num_targets_in_cluster_with_label
            table_d = num_targets - num_targets_with_label - table_b

            oddsratio, pval = stats.fisher_exact([[table_a, table_b], [table_c, table_d]])

            num_tests = num_labels_in_cluster
            corrected_pval = num_tests * pval

            label_freq = num_targets_with_label/num_targets
            label_freq_in_cluster = num_targets_in_cluster_with_label/num_targets_in_cluster

            rows.append(
                Result(
                    cluster_id=cluster_id,
                    label_id=label_id,
                    cluster_size=num_targets_in_cluster,
                    label_enrichment=(label_freq_in_cluster/label_freq),
                    corrected_pval = corrected_pval,
                )
            )

    results = (
        pd.DataFrame(rows)
        .sort_values(by=['cluster_id', 'label_id', 'corrected_pval'], ascending=True)
    )
    return results