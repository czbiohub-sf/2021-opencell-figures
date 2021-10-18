import seaborn as sns

from ..external import sankey
from .definitions import label_orders, consensus_label_colors


def plot_sankey(
    df, 
    res='high', 
    left_category_name='OpenCell', 
    right_category_name='', 
    min_count=10, 
    use_dynamic_colormap=True
):
    '''
    Compare two sets of labels using a sankey diagram
    df : dataframe of annotations (one row per target)
    min_count : the minimum number of times a label set must occur to be included in the diagram
        (this threshold is applied to both the left_label and right_label columns)
    '''

    # separator to use for concatenating label sets
    sep = ', '

    # the column names containing the categories to show on the left and right of the sankey diagram
    left_label = 'consensus_label_oc'
    right_label = (
        'consensus_label_hpa' if 'consensus_label_hpa' in df.columns else 'consensus_label_yeast'
    )

    def concat_labels(labels):
        return sep.join(sorted(labels))

    # construct readable sankey labels by concatenating the label sets
    df = df.copy()
    df[left_label] = df[left_label].apply(concat_labels)
    df[right_label] = df[right_label].apply(concat_labels)

    sankey_label_order = [concat_labels(labels) for labels in label_orders[res][::-1]]

    # find the most common label sets among both HPA and OC labels
    most_common_label_sets = []
    counts = df[left_label].value_counts()
    most_common_label_sets.extend(counts[counts > min_count].index.tolist())

    counts = df[right_label].value_counts()
    most_common_label_sets.extend(counts[counts > min_count].index.tolist())

    # we only consider targets whose OC *and* HPA label sets are *both* among the most common
    mask = (
        df[left_label].isin(most_common_label_sets) & df[right_label].isin(most_common_label_sets)
    )

    df = df[mask]
    left_label_order = [label for label in sankey_label_order if label in df[left_label].values]
    right_label_order = [label for label in sankey_label_order if label in df[right_label].values]

    all_label_sets = set(left_label_order).union(right_label_order)
    if use_dynamic_colormap:
        colormap = create_dynamic_colormap(all_label_sets, label_set_order=sankey_label_order)
    else:
        colormap = create_hardcoded_colormap(all_label_sets, res=res, sep=sep)

    sankey.sankey(
        df[left_label], 
        df[right_label], 
        leftLabels=left_label_order,
        rightLabels=right_label_order,
        left_category_name=left_category_name,
        right_category_name=right_category_name,
        colorDict=colormap, 
        aspect=10, 
        fontsize=11
    )


def create_hardcoded_colormap(label_sets, res, sep):
    '''
    Calculate a colormap for each label set using a hardcoded colormap for the consensus labels
    (imported as consensus_label_colors)
    For label sets with two labels (which is most of them), we average the two labels' colors
    For label sets with three labels, we do nothing (use gray)
    '''
    colormap = {}
    for label_set in label_sets:

        labels = label_set.split(sep)
        if len(labels) == 1:
            color = consensus_label_colors[res][labels[0]]

        # interpolate between the two colors
        elif len(labels) == 2:
            color = sns.blend_palette(
                [consensus_label_colors[res][labels[0]], consensus_label_colors[res][labels[1]]],
                n_colors=3
            )[1]

        # too many labels to use a color, so use gray
        else:
            color = '#878584'

        colormap[label_set] = color
    return colormap


def create_dynamic_colormap(label_sets, label_set_order):
    '''
    Calculate a colormap for the label_sets from a hard-coded HSV-like colormap,
    usin the order defined by label_set_order
    (label_sets is assumed to be a subset of the label sets in label_set_order)

    The motivation for this is to make pretty sankey diagrams 
    in which adjacent colors are adjacent in the HSV colormap.

    Note that there are too many label_sets in label_set_order to define one global colormap,
    and besides, when label_sets is a small subset of label_set_order, 
    a global colormap would result in large jumps in hue between adjacent labels in the sankey.

    Note that this works best if len(label_sets) is less than len(hsv_colors); 
    if it's not, then the colormap is cyclic and less effective.
    '''
    hsv_colors = [
        '#f44336', '#e91e63', '#9c27b0', '#673ab7', '#3f51b5',
        '#2196f3', '#03a9f4', '#00bcd4', '#009688', '#4caf50',
        '#8bc34a', '#cddc39', '#ffeb3b', '#ffc107', '#ff9800',
        '#ff5722', '#795548', '#009e9e', '#607d8b', 
    ]

    # colormap = {
    #     label_set: hsv_colors[ind % len(hsv_colors)]
    #     for ind, label_set in enumerate(label_set_order)
    # }
    # return colormap

    colormap = {}
    counter = 0
    for label_set in label_set_order:
        if label_set in label_sets:
            colormap[label_set] = hsv_colors[counter % len(hsv_colors)]
            counter += 1
    
    return colormap
