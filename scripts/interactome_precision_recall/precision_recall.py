import pandas as pd
import numpy as np

from itertools import repeat
from multiprocessing import Pool
from pyseus import validation_analysis as vali


def precision_recall_curve(all_hits, target_col, prey_col, corum, helas, thresholds, metric=None):
    """
    calculates colocalization precision and corum coverage recall to create precision-recall
    curve data by sliding across interaction calling thresholds

        all_hits – DataFrame, output of get_all_interactors or drop_unchosen, just_hits=False
        corum – DataFrame, CORUM complex dataframe
        helas – DataFrame, protein group localization data from fractions MS, despite the naming,
        it can either be HELA or HEK dataset
        thresholds – list, a list of floats that serve as thresholds for interaction calling
        custom – boolean, whether to use simple p-value thresholding or FDR curve thresholding

    """

    # if a custom metric is given, use the column, otherwise use standard pval-enrichment cols
    if metric:
        all_hits = all_hits[[target_col, prey_col, metric]].copy()
    else:
        all_hits = all_hits[[target_col, prey_col , 'pvals', 'enrichment']].copy()
    all_hits['target'] = all_hits['target'].astype(str)
    all_hits['prey'] = all_hits['prey'].astype(str)

    recalls = []
    precisions = []
    
    # For recall-precisions, set default threshold if threshold is not defined
    if not thresholds:
        if metric:
            thresholds = [1, 2, 3, 4.5, 6, 9, 12, 15, 20, 25, 30]
        else:
            thresholds = [0, 1, 2, 3, 4, 5, 6, 7]

    multi_args = zip(repeat(all_hits), repeat(corum), repeat(helas), repeat(metric), thresholds)

    # multi processing
    p = Pool()
    precision_recalls = p.starmap(precision_recall, multi_args)

    p.close()
    p.join()

    unzip_prs = list(zip(*precision_recalls))
    precisions, recalls = unzip_prs[0], unzip_prs[1]
    pr_table = pd.DataFrame()
    pr_table['precision'] = precisions
    pr_table['recall'] = recalls

    return pr_table


def precision_recall(all_hits, corum, helas, metric, threshold, interaction_called=False):
    """
    function for calulating precision recall in single threshold, used for parallel processing
    """

    all_hits = all_hits.copy()

    if interaction_called:
        analysis = vali.Validation(
            hit_table=None,
            target_col='target', 
            prey_col='prey',
            corum=corum, 
            localization_table=helas, 
            interaction_table=all_hits
        )
    else:
        analysis = vali.Validation(
            hit_table=all_hits, 
            target_col='target',
            prey_col='prey',
            corum=corum, 
            localization_table=helas
        )
        if metric:
            analysis.pval_threshold(threshold, metric=metric)
        else:
            analysis.static_fdr(curvature=3, offset=threshold)

    # calculate recall
    analysis.corum_interaction_coverage()
    recall = analysis.recall

    # calculate precision
    analysis.colocalization_precision()
    precision = analysis.precision

    return [precision, recall]


