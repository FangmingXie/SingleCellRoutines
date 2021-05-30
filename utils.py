"""
"""
from __init__scr import *

import numpy as np
import pandas as pd
import logging
import subprocess as sp
import os
import collections
from concurrent.futures import TimeoutError

def create_logger(name='log'):
    """
    args: logger name
    return: a logger object
    """
    logging.basicConfig(
        format='%(asctime)s %(message)s', 
        level=logging.INFO)
    return logging.getLogger(name)

def task_done(future):
    """parallelize utils
    """
    try:
        result = future.result()  # blocks until results are ready
    except TimeoutError as error:
        print("Function took longer than %d seconds" % error.args[1])
    except Exception as error:
        print("Function raised %s" % error)
        print(error.traceback)  # traceback of the function

def compress(file, suffix='gz'):
    """
    """
    # bgzip compression
    try:
        sp.run("bgzip -f {}".format(file), shell=True)
        if suffix != 'gz':
            sp.run("mv {}.gz {}.{}".format(file, file, suffix), shell=True)
    except:
        sp.call("bgzip -f {}".format(file), shell=True)
        if suffix != 'gz':
            sp.call("mv {}.gz {}.{}".format(file, file, suffix), shell=True)
    return

def import_single_textcol(fname, header=None, col=0):
    return pd.read_csv(fname, header=header, sep='\t')[col].values

def export_single_textcol(fname, array):
    with open(fname, 'w') as f:
        f.write('\n'.join(array)+'\n')

def get_mCH_contexts():
    """
    """
    contexts = []
    for base1 in ['A','C','T']:
        for base2 in ['A','C','G','T']:
            contexts.append('C' + base1 + base2)
    return contexts+['CAN', 'CTN', 'CCN']

def get_expanded_context(context):
    """
    """
    if context == "CH":
        # 15 contexts 
        contexts = get_mCH_contexts()
    elif context == "CG":
        contexts = ["CGA","CGC","CGG","CGT","CGN"]
    elif context == "CA":
        contexts = ["CAA","CAC","CAG","CAT","CAN"]
    elif context == "CT":
        contexts = ["CTA","CTC","CTG","CTT","CTN"]
    else: 
        contexts = [context]
    return contexts

def tabix_summary(records, context="CH", cap=0):
    """
    """
    mc = 0
    c = 0
    contexts = get_expanded_context(context)
    if cap > 0:
        for record in records:
            if record[3] in contexts:
                if int(record[5]) <= cap:
                    mc += int(record[4])
                    c += int(record[5])
    else:
        for record in records:
            if record[3] in contexts:
                mc += int(record[4])
                c += int(record[5])
    return mc, c

def read_allc(fname, pindex=True, compression='gzip', remove_chr=True, **kwargs):
    """
    """
    if pindex:
        df = pd.read_csv(fname, 
            sep="\t", 
            compression=compression,
            header=None, 
            index_col=['chr', 'pos'],
            dtype={'chr': str, 'pos': np.int, 'mc': np.int, 'c': np.int, 'methylated': np.int},
            names=['chr','pos','strand','context','mc','c','methylated'], **kwargs)
    else:
        df = pd.read_csv(fname, 
            sep="\t", 
            compression=compression,
            header=None, 
            # index_col=['chr', 'pos'],
            dtype={'chr': str, 'pos': np.int, 'mc': np.int, 'c': np.int, 'methylated': np.int},
            names=['chr','pos','strand','context','mc','c','methylated'], **kwargs)
        
    # remove chr
    if remove_chr:
        if df.iloc[0,0].startswith('chr'):
            df['chr'] = df['chr'].apply(lambda x: x[3:]) 
    return df

def get_chrom_lengths(genome_size_file):
    """
    """
    chrom_sizes = pd.read_csv(genome_size_file, sep="\t", header=None) 
    chrom_sizes = chrom_sizes[~chrom_sizes[0].str.contains('_')]
    # remove leading 'chr'
    chrom_sizes[0] = chrom_sizes[0].apply(lambda x: x[len('chr'):])
    chrom_sizes = chrom_sizes.sort_values(0).set_index(0).squeeze()
    return chrom_sizes

def anova_analysis(df, groups, per_group_tss=False):
    """ANOVA analysis
    df is a dataframe (n_obs, n_features)
    group_assigments is an array (n_obs,) that matches df
    """
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    N = len(df) 
    K = len(np.unique(groups))

    mean = df.mean()
    mean.name = 'mu'

    tss = np.power(df-mean, 2).sum()
    tss.name = 'tss'

    tmp = df.copy()
    tmp['_group'] = groups

    # tss_in
    tss_in = 0
    for grp, dfsub in tmp.groupby('_group'):
        diff = dfsub.drop('_group', axis=1) - dfsub.drop('_group', axis=1).mean()
        tss_in += np.power(diff, 2).sum()
    tss_in.name = 'tss_in'

    # res
    tss_df = tss.to_frame().join(tss_in)
    tss_df['tss_out'] = tss_df['tss'] - tss_df['tss_in']
    tss_df['eta2'] = tss_df['tss_out']/tss_df['tss']
    tss_df['s2'] = tss_df['tss']/len(df)
    tss_df = tss_df.join(mean)

    # test 
    tss_df['f'] = (tss_df['tss_out']/(K-1))/(tss_df['tss_in']/(N-K))
    tss_df['p_f'] = stats.f.sf(tss_df['f'], N, K)
    tss_df['p_f'] = tss_df['p_f'].fillna(1)
    rej, fdr, _, _ = multipletests(tss_df['p_f'].values, method='fdr_bh')
    tss_df['fdr_f'] = fdr 
    tss_df['-log10fdr_f'] = -np.log10(fdr)

    # show individual contributions
    if per_group_tss:
        group_mean = tmp.groupby('_group').mean()
        group_size = tmp.groupby('_group').size()
        group_tss = np.power((group_mean - mean), 2).multiply(group_size, axis=0)
        group_eta2_frac = group_tss/tss_df['tss_out']

        group_eta2_frac.index = ['eta2_frac_'+group for group in group_tss.index]
        tss_df = tss_df.join(group_eta2_frac.T)

    return tss_df 

def get_order_from_hierarchy(mat, **kwargs):
    """(n_obs, n_features)
    """
    import scipy.cluster.hierarchy as sch

    Z = sch.linkage(mat, **kwargs)
    d = sch.dendrogram(Z, no_plot=True)
    order = d['leaves']

    return order

def get_fdr(p, method='fdr_bh', **kwargs):
    """
    """
    from statsmodels.stats.multitest import multipletests

    _, fdr, _, _ = multipletests(p, method=method, **kwargs)
    return fdr

def savefig(fig, path, dpi=300):
    """
    """
    fig.savefig(path, bbox_inches='tight', dpi=300)
    return 