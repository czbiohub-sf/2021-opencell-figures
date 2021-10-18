import urllib.parse
import urllib.request
import sys
import multiprocessing
import os
import re
import pickle
import pandas as pd
import numpy as np
from itertools import repeat
from multiprocessing import Pool


class RawTables:
    """
    Raw Tables class contains DataFrame objects, functions, and metadata that cover
    multiple pre-processing steps to create a final processed imputed table. 
    
    """
    
    # initiate raw table by importing from data directory
    def __init__(self, root, analysis, intensity_type, pg_file='proteinGroups.txt'):
        
        self.raw_table = pd.read_csv(root+pg_file,
            sep='\t', header=0, low_memory=False)
        
        # root directory
        self.root = root
        # analysis string
        self.analysis = analysis
        # Specification of which intensity (raw or LFQ) to use
        self.intensity_type = intensity_type
    

    def save(self, option_str=''):
        """
        save class to a designated directory
        """
        analysis_dir = self.root + self.analysis
        if len(option_str) > 0:
            option_str = '_' + option_str
        file_dir = analysis_dir + "/preprocessed_tables" + option_str + '.pkl' 
        if not os.path.isdir(analysis_dir):
            print(analysis_dir)
            print('Directory does not exist! Creating new directory')
            os.mkdir(analysis_dir)

        print("Saving to: " + file_dir)
        with open(file_dir, 'wb') as file_:
            pickle.dump(self, file_, -1)

      
    def filter_table(self, verbose=True):
        """filter rows that do not meet the QC (contaminants, reverse seq, only identified by site)
        Also filter non-intensity columns that will not be used for further processing"""
        
        ms_table = self.raw_table.copy()

        pre_filter = ms_table.shape[0]

        # remove rows with potential contaminants
        ms_table = ms_table[ms_table['Potential contaminant'].isna()]

        # remove rows only identified by site
        ms_table = ms_table[ms_table['Only identified by site'].isna()]

        # remove rows that are reverse seq
        ms_table = ms_table[ms_table['Reverse'].isna()]

        filtered = pre_filter - ms_table.shape[0]
        if verbose:
            print("Filtered " + str(filtered) + ' of '
                  + str(pre_filter) + ' rows. Now '
                  + str(ms_table.shape[0]) + ' rows.')
        
        # select necessary columns
        all_cols = list(ms_table)
        info_cols = ['Protein IDs', 'Majority protein IDs', 'Protein names',
        'Gene names']
        intensity_cols = select_intensity_cols(all_cols, self.intensity_type)

        ms_table = ms_table[info_cols + intensity_cols]


        self.filtered_table = ms_table
    
    def transform_intensities(self, func=np.log2):
        """transform intensity values in the dataframe to a given function"""
        
        try:
            filtered = self.filtered_table.copy()
        except AttributeError:
            print(
                "Raw table has not been filtered yet, use filter_table() method"\
                "before transforming intensities")
            return

        filtered = self.filtered_table.copy()
        intensity_cols = select_intensity_cols(list(self.filtered_table),
            intensity_type=self.intensity_type)

        # for each intensity column, transform the values
        for int_col in intensity_cols:
            # if transformation is log2, convert 0s to nans
            # (faster in one apply step than 2)
            if func == np.log2:
                filtered[int_col] = filtered[int_col].apply(lambda x: np.nan
                    if x == 0 else func(x))
            else:
                filtered[int_col] = filtered[int_col].apply(func)
                # Replace neg inf values is np.nan
                filtered[int_col] = filtered[int_col].apply(
                    lambda x: np.nan if np.isneginf(x) else x)
        
        self.transformed_table = filtered
    
    def group_replicates(self, intensity_re=r'_\d+$', reg_exp=r'(.*_.*)_\d+$'):
        """Group the replicates of intensities into replicate groups"""
        
        reg_exp = self.intensity_type + ' ' + reg_exp
        try: 
            self.transformed_table
            transformed = self.transformed_table.copy()
        except AttributeError:
            print(
            'Intensity values have not been transformed yet from '\
            'filtered table,\nwe recommend using transform_intensities() '\
            'method before grouping replicates.\n')

            try: 
                self.filtered_table
                print("Using filtered_table to group replicates.")
                transformed = self.filtered_table.copy()
            except AttributeError:
                print('Please filter raw table first using filter_table()\
                    method.')
                return


       # get col names
        col_names = list(transformed)

        # using a dictionary, group col names into replicate groups
        group_dict = {}
        for col in col_names:
            # search REs of the replicate ID, and get the group names

            # search if the col is for intensity values
            intensity_search = re.search(intensity_re, col.lower(),
                flags=re.IGNORECASE)

            # if so, get the group name and add to the group dict
            # use groups from re.search to customize group names
            if intensity_search:
                group_search = re.search(reg_exp, col, flags=re.IGNORECASE)
                group_name = ''

                for re_group in group_search.groups():
                    group_name += re_group
                group_dict[col] = group_name

            # if not, group into 'Info'
            else:
                group_dict[col] = 'Info'


        # pd function to add the replicate group to the columns
        grouped = pd.concat(dict((*transformed.groupby(group_dict, 1),)), axis=1)

        grouped.columns = grouped.columns.rename("Baits", level=0)
        grouped.columns = grouped.columns.rename("Replicates", level=1) 

        self.grouped_table = grouped
    
    
    def remove_invalid_rows(self):
        """Remove rows that do not have at least one group that has values
        in all triplicates"""

        try:
            grouped = self.grouped_table.reset_index(drop=True).copy()
        except AttributeError:
            print("Replicates need to be grouped before this method."\
                "Please use group_replicates() to group replicates under same sample")
            return

        # reset index
        grouped = self.grouped_table.reset_index(drop=True).copy()
        unfiltered = self.grouped_table.shape[0]

        # Get a list of all groups in the df
        group_list = list(set([col[0] for col in list(grouped) if col[0] != 'Info']))

        # booleans for if there is a valid value
        filtered = grouped[group_list].apply(np.isnan)
        # loop through each group, and filter rows that have valid values
        for group in group_list:
            # filter all rows that qualify as all triplicates having values
            filtered = filtered[filtered[group].any(axis=1)]

        # a list containing all the rows to delete
        del_list = list(filtered.index)

        # create a new df, dropping rows with invalid data
        filtered_df = grouped.drop(del_list)
        filtered_df.reset_index(drop=True, inplace=True)
        filtered = filtered_df.shape[0]

        print("Removed invalid rows. " + str(filtered) + " from "
            + str(unfiltered) + " rows remaining.")

        self.preimpute_table = filtered_df
    
    def bait_impute(self, distance=1.8, width=0.3, local=True):
        """
        bait-imputation for sets of data without enough samples.
        This fx imputes a value from a normal distribution of the left-tail
        of a bait’s capture distribution for the undetected preys using
        multi-processing.
            distance: float, distance in standard deviation from the
            mean of the sample distribution upon which to impute. Default = 0
            width: float, width of the distribution to impute in standard deviations. Default = 0.3
        """
        
        try:
            imputed = self.preimpute_table.copy()
        except AttributeError:
            print("group_replicates() and remove_invalid_rows() need to be run"\
                "before imputation")
            return
        
        self.bait_impute_params = {'distance': distance, 'width': width}

        # Retrieve all col names that are not classified as Info
        bait_names = [col[0] for col in list(imputed) if col[0] != 'Info']
        baits = list(set(bait_names))
        bait_series = [imputed[bait].copy() for bait in baits]
        if local:
            global_mean = 0
            global_stdev = 0

        else:
            # if not using columnwise imputation, calculate global mean and stdev
            all_intensities = imputed[[col for col in list(imputed) if col[0] != 'Info']].copy()
            global_mean = all_intensities.droplevel('Baits', axis=1).stack().mean()
            global_stdev = all_intensities.droplevel('Baits', axis=1).stack().std()



        bait_params = zip(
            bait_series, repeat(distance), repeat(width), repeat(local), repeat(global_mean), repeat(global_stdev))

        # Use multiprocessing pool to parallel impute
        p = Pool()
        impute_list = p.starmap(pool_impute, bait_params)
        p.close()
        p.join()

        for i, bait in enumerate(baits):
            imputed[bait] = impute_list[i]

        self.bait_imputed_table = imputed
    
    def prey_impute(self, distance=0, width=0.3, thresh=100):
        """
        default mode of imputation. For protein groups with less than threshold number
        of sample number, impute a value from a normal distribution of the prey’s capture
        distribution using multi-processing. Note- most protein groups do not need imputation
        with 12-plate MBR

            distance: float, distance in standard deviation from the mean of the
                sample distribution upon which to impute. Default = 0
            width: float, width of the distribution to impute in standard deviations.
                Default = 0.3
            threshold: int, max number of samples required for imputation
        """
        
        try:
            imputed = self.preimpute_table.copy()
        except AttributeError:
            print("group_replicates() and remove_invalid_rows() need to be run"\
                "before imputation")
            return
        
        imputed = self.preimpute_table.copy()
        imputed.drop(columns='Info', inplace=True)
        imputed = imputed.T
        self.prey_impute_params = {'distance': distance, 'width': width,
            'thresh': thresh}

        # Retrieve all col names that are not classified as Info
        baits = list(imputed)
        bait_series = [imputed[bait].copy() for bait in baits]
        bait_params = zip(
            bait_series, repeat(distance), repeat(width), repeat(thresh))

        # Use multiprocessing pool to parallel impute
        p = Pool()
        impute_list = p.starmap(pool_impute_prey, bait_params)
        p.close()
        p.join()

        for i, bait in enumerate(baits):
            imputed[bait] = impute_list[i]

        imputed = imputed.T

        info_cols = [x for x in list(self.preimpute_table) if x[0] == 'Info']
        for col in info_cols:
            imputed[col] = self.preimpute_table[col]      

        self.prey_imputed_table = imputed  
    

    def generate_export_bait_matrix(self):
        """
        Generates and creates a Boolean bait matrix that will be used for control
        exclusion in p-val and enrichment analysis. 
        """
        grouped = self.grouped_table.copy()
        baits = list(set(grouped.columns.get_level_values('Baits').to_list()))
        baits.remove('Info')
        baits.sort()
        bait_df = pd.DataFrame()
        bait_df['Baits'] = baits
        bait_df.reset_index(drop=True, inplace=True)
        self.bait_list = bait_df.copy()
        bait_df2 = bait_df.copy()
        bait_df2['plot'] = True
        bait_df2.to_csv(self.root + self.analysis + '/plotting_exclusion_list.csv',
            index=False)

        # Create a boolean table
        for bait in baits:
            bait_df[bait] = True
        self.bait_matrix = bait_df.copy()
        self.bait_matrix.to_csv(self.root + self.analysis + '/analysis_exclusion_matrix.csv',
            index=False)


def czb_initial_processing(root, analysis, pg_file='proteinGroups.txt',
    intensity_type='LFQ intensity', bait_impute=True, distance=1.8, width=0.3,
    thresh=100, local=True):
    
    """
    wrapper script for all the pre-processing up to imputation using
    PyseusRawTables Class. Saves and returns the PyseusRawTables in the
    designated analysis directory
    """
    # make directory for analysis folder
    analysis_dir = root + analysis

    if not os.path.isdir(analysis_dir):
        os.mkdir(analysis_dir)
    
    # Run all the processing methods
    pyseus_tables = RawTables(root=root,
        analysis=analysis, intensity_type=intensity_type, pg_file=pg_file)
    pyseus_tables.filter_table()
    pyseus_tables.transform_intensities(func=np.log2)
    pyseus_tables.group_replicates(intensity_re=r'_\d+$', reg_exp=r'(.*_.*)_\d+$')
    pyseus_tables.remove_invalid_rows()
    if bait_impute:
        pyseus_tables.bait_impute(distance=distance, width=width, local=local)
    else:
        pyseus_tables.prey_impute(distance=distance, width=width, thresh=thresh)
    pyseus_tables.generate_export_bait_matrix()
    pyseus_tables.save()
    return pyseus_tables


def load_raw_tables(file_dir):
    """
    use pickle to load RawTables class
    """
    return pickle.load(open(file_dir, 'rb', -1))

        
def select_intensity_cols(orig_cols, intensity_type):
    """from table column names, return a list of only intensity cols
    rtype: intensity_cols list """
    # new list of intensity cols
    intensity_cols = []

    # create a regular expression that can distinguish between
    # intensity and LFQ intensity
    re_intensity = '^' + intensity_type.lower()

    # for loop to include all the intensity col names
    intensity_type = intensity_type.lower()
    for col in orig_cols:
        col_l = col.lower()

        # check if col name has intensity str
        if re.search(re_intensity, col_l):
            intensity_cols.append(col)

    return intensity_cols


def pool_impute(bait_group, distance=1.8, width=0.3, local=True, global_mean=0, global_stdev=0):
    """target for multiprocessing pool from multi_impute_nans"""
    all_vals = bait_group.stack()
    mean = all_vals.mean()
    stdev = all_vals.std()
    
    if local:
        # get imputation distribution mean and stdev
        imp_mean = mean - distance * stdev
        imp_stdev = stdev * width
    else:
        # use global mean and stdev
        imp_mean = global_mean
        imp_stdev = global_stdev


    # copy a df of the group to impute values
    bait_df = bait_group.copy()

    # loop through each column in the group
    for col in list(bait_df):
        bait_df[col] = bait_df[col].apply(random_imputation_val,
            args=(imp_mean, imp_stdev))
    return bait_df


def pool_impute_prey(bait_group, distance=0, width=0.3, thresh=100):
    """target for multiprocessing pool from multi_impute_nans"""

    if bait_group.count() > thresh:
        return bait_group


    mean = bait_group.mean()
    stdev = bait_group.std()

    # get imputation distribution mean and stdev
    imp_mean = mean - distance * stdev
    imp_stdev = stdev * width

    # copy a df of the group to impute values
    bait_df = bait_group.copy()


    bait_df = bait_df.apply(random_imputation_val,
            args=(imp_mean, imp_stdev))
    return bait_df


def random_imputation_val(x, mean, std):
    """from a normal distribution take a random sample if input is
    np.nan. For real values, round to 4th decimal digit.
    Floats with longer digits will be 'barcoded' by further digits

    rtype: float"""

    if np.isnan(x):
        return np.random.normal(mean, std, 1)[0]
    else:
        return np.round(x, 4)

def sample_rename(col_names, RE, replacement_RE, repl_search=False):
    """
    method to change column names for previewing in notebook
    """

    # start a new col list
    new_cols = []

    # Loop through cols and make quaifying subs
    for col in col_names:
        for i in np.arange(len(RE)):
            if re.search(RE[i], col, flags=re.IGNORECASE):
                replacement = replacement_RE[i]
                if (repl_search) & (len(replacement) > 1):
                    rep_search = re.search(replacement, col,
                                flags=re.IGNORECASE)
                    replacement = ''
                    for group in rep_search.groups():
                        replacement += group

                col = re.sub(RE[i], replacement, col, flags=re.IGNORECASE)
        new_cols.append(col)
    return new_cols


def rename_columns(df, RE, replacement_RE, repl_search=False):
    """
    change intensity column names to a readable format. More specifically,
    search a column name from an input RE and substitute matches with another
    input substitute strings or REs.
        col_names: list, a list of column names from raw_df
        RE: list, a list of regular expressions to search in column names
        replacement_RE: list, a list of strs/REs that substitute the original expression
        repl_search: boolean, if True, elements in replacement_RE are treated as regular
            expressions used in search, and all specified groups are used in substitution

    """
    df = df.copy()
    col_names = list(df)

    # start a new col list
    new_cols = []

    # Loop through cols and make quaifying subs
    for col in col_names:
        for i in np.arange(len(RE)):
            if re.search(RE[i], col, flags=re.IGNORECASE):
                replacement = replacement_RE[i]
                if (repl_search) & (len(replacement) > 1):
                    rep_search = re.search(replacement, col,
                                flags=re.IGNORECASE)
                    replacement = ''
                    for group in rep_search.groups():
                        replacement += group

                col = re.sub(RE[i], replacement, col, flags=re.IGNORECASE)
        new_cols.append(col)
    
    rename = {i: j for i, j in zip(col_names, new_cols)}

    renamed = df.rename(columns=rename)

    return renamed


def median_replicates(imputed_df, mean=False, save_info=True, col_str=''):
    """For each bait group, calculate the median of the replicates
    and returns a df of median values

    rtype: median_df pd dataframe"""

    imputed_df = imputed_df.copy()
    # retrieve bait names
    bait_names = [col[0] for col in list(imputed_df) if col[0] != 'Info']
    bait_names = list(set(bait_names))
    # initiate a new df for medians
    median_df = pd.DataFrame()

    # for each bait calculate medain across replicates and add
    # to the new df
    for bait in bait_names:
        if mean:
            bait_median = imputed_df[bait].mean(axis=1)
        else:
            bait_median = imputed_df[bait].median(axis=1)
        new_col_name = col_str + bait
        median_df[new_col_name] = bait_median

    if save_info:
        # get info columns into the new df
        info = imputed_df['Info']
        info_cols = list(info)

        for col in info_cols:
            median_df[col] = info[col]

    return median_df