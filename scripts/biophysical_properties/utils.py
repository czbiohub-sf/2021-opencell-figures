import pandas as pd
import numpy as np
import requests
import dask
import dask.diagnostics
import metapredict as mp

from Bio.SeqUtils.ProtParam import ProteinAnalysis


def get_sequence_from_uniprot(uniprot_id):
    '''
    Retrieve a protein sequence from the uniprot API
    '''
    url = f'https://www.uniprot.org/uniprot/{uniprot_id}.fasta'
    response = requests.get(url)
    try:
        response.raise_for_status()
    except:
        return None
    return response.text


def get_sequences_from_uniprot(uniprot_ids):
    '''
    Retrieve the amino acid sequences for a list of uniprot_ids
    Returns a dataframe with columns uniprot_id, sequence, fasta
    '''
    tasks = [dask.delayed(get_sequence_from_uniprot)(uniprot_id) for uniprot_id in uniprot_ids]
    with dask.diagnostics.ProgressBar():
        fastas = dask.compute(*tasks)

    seqs = pd.DataFrame(data={'uniprot_id': uniprot_ids, 'fasta': fastas})

    seqs = seqs.loc[~seqs.fasta.isna()]
    seqs['sequence'] = seqs.fasta.apply(lambda s: ''.join(s.split('\n')[1:]))

    missing_uniprot_ids = set(uniprot_ids).difference(seqs.uniprot_id)
    if len(missing_uniprot_ids):
        print('No results for %s uniprot_ids' % len(missing_uniprot_ids))

    return seqs, missing_uniprot_ids


def _get_iupred_scores(uniprot_id):
    '''
    Get the IUPred2A disorder scores for a single uniprot_id using the iupred API
    '''
    try:
        response = requests.get('http://iupred2a.elte.hu/iupred2a/%s' % uniprot_id)
        df = pd.DataFrame(data=[line.split('\t') for line in response.text.split('\n')[5:]])
        df = df.dropna(how='any', axis=0)
        df.columns = ['position', 'aa', 'score']
        df['uniprot_id'] = uniprot_id
    except:
        df = pd.DataFrame(data=[{'uniprot_id': uniprot_id}])
    return df


def get_iupred_scores(uniprot_ids):
    '''
    Get the IUPred2A disorder scores for each uniprot_id in `uniprot_ids`
    '''
    tasks = [dask.delayed(_get_iupred_scores)(uniprot_id) for uniprot_id in uniprot_ids]
    with dask.diagnostics.ProgressBar():
        dfs = dask.compute(*tasks)

    all_dfs = pd.concat(dfs, axis=0)
    all_dfs.dropna(axis=0, how='any', inplace=True)

    # missing uniprot_ids
    missing_ids = set(uniprot_ids).difference(all_dfs.uniprot_id)
    print('No scores retrieved for %s uniprot_ids; trying again' % len(missing_ids))

    # try again to retrieve the missing ids
    tasks = [dask.delayed(_get_iupred_scores)(uniprot_id) for uniprot_id in missing_ids]
    with dask.diagnostics.ProgressBar():
        dfs = dask.compute(*tasks)

    all_dfs = pd.concat((all_dfs, *dfs), axis=0)
    all_dfs.dropna(axis=0, how='any', inplace=True)

    all_dfs['score'] = all_dfs.score.astype(float)
    all_dfs['position'] = all_dfs.position.astype(int)

    all_dfs = (
        all_dfs
        .groupby(['uniprot_id', 'position'])
        .first()
        .reset_index()
        .sort_values(by=['uniprot_id', 'position'])
    )
    return all_dfs


def calc_metapredict_scores(sequences):
    '''
    Calculate metapredict scores for a given set of sequences

    sequences : dataframe with columns 'uniprot_id', 'sequence'
    Returns : a dataframe of metapredict scores with one row per sequence position
        and columns uniprot_id, position, score
    '''
    def _calc_scores(row):
        try:
            scores = mp.predict_disorder(row.sequence)
        except:
            scores = [np.nan]

        df = pd.DataFrame(data={'score': scores, 'position': np.arange(len(scores))})
        df['uniprot_id'] = row.uniprot_id
        return df

    tasks = [dask.delayed(_calc_scores)(row) for ind, row in sequences.iterrows()]
    with dask.diagnostics.ProgressBar():
        dfs = dask.compute(*tasks)
    
    scores = pd.concat(dfs, axis=0)
    return scores


def calc_windowed_disorder_scores(scores, window_sizes):
    '''
    Calculate the maximum mean disorder score in sliding windows of various sizes,
    using either the IUPred2A or the metapredict scores

    scores : dataframe of scores returned by calc_metapredict_scores
    windows : list of window sizes (in number of aas)
    '''
    uniprot_ids = scores.uniprot_id.unique()

    grouped = scores.groupby('uniprot_id')
    def _calc_windowed_scores(uniprot_id, window_sizes):
        score = grouped.score.get_group(uniprot_id)
        overall_mean = score.mean()
        max_windowed_means = [
            score.rolling(window_size).mean().max() for window_size in window_sizes
        ]
        return [uniprot_id, overall_mean] + max_windowed_means

    tasks = [
        dask.delayed(_calc_windowed_scores)(uniprot_id, window_sizes)
        for uniprot_id in uniprot_ids[:]
    ]
    with dask.diagnostics.ProgressBar():
        rows = dask.compute(*tasks)

    columns = (
        ['uniprot_id', 'mean_score'] + 
        [f'max_window_{window_size}aa' for window_size in window_sizes]
    )
    windowed_scores = pd.DataFrame(data=rows, columns=columns)
    return windowed_scores


def _calc_biophysical_properties(sequence, uniprot_id=None):
    '''
    Calculate basic biophysical properties for an amino acid sequence
    using biopython's ProteinAnalysis class
    '''
    props = dict(uniprot_id=uniprot_id)
    protein_analysis = ProteinAnalysis(sequence)

    # properties that are direct attributes of the ProtParam.ProteinAnalysis class
    attrs = ['molecular_weight', 'aromaticity', 'instability_index', 'isoelectric_point', 'gravy']
    for attr in attrs:
        props[attr] = getattr(protein_analysis, attr)()
    
    # properties that require an argument or averaging
    props['charge_at_neutral_ph'] = protein_analysis.charge_at_pH(7.0)
    props['mean_flexibility'] = np.mean(protein_analysis.flexibility())

    # fraction of amino acids that tend to be in helix, turns, sheets
    helix, turn, sheet = protein_analysis.secondary_structure_fraction()
    props['fraction_helix'] = helix
    props['fraction_turn'] = turn
    props['fraction_sheet'] = sheet
    
    return props


def calc_biophysical_properties(sequences):
    '''
    sequences : dataframe of protein sequences with columns uniprot_id, sequence
    '''

    def _calc_props(sequence, uniprot_id):
        try:
            props = _calc_biophysical_properties(sequence, uniprot_id)
        except:
            props = dict(uniprot_id=uniprot_id)
        return props

    tasks = [
        dask.delayed(_calc_props)(row.sequence, row.uniprot_id)
        for ind, row in sequences.iterrows()
    ]
    with dask.diagnostics.ProgressBar():
        rows = dask.compute(*tasks)
    
    props = pd.DataFrame(data=rows)
    return props
