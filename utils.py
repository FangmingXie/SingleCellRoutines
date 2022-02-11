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

def diag_matrix(X, rows=np.array([]), cols=np.array([]), threshold=None):
    """Diagonalize a matrix as much as possible
    threshold controls the level of diagnalization
    a smaller threshold enforces more number of strict diagnal values,
    while encourages less number of free columns (quasi-diagnal)
    """
    di, dj = X.shape
    transposed = 0
    
    # enforce nrows <= ncols
    if di > dj:
        di, dj = dj, di
        X = X.T.copy()
        rows, cols = cols.copy(), rows.copy()
        transposed = 1
        
    # start (di <= dj)
    new_X = X.copy()
    new_rows = rows.copy() 
    new_cols = cols.copy() 
    if new_rows.size == 0:
        new_rows = np.arange(di)
    if new_cols.size == 0:
        new_cols = np.arange(dj)
        
    # bring the greatest values in the lower right matrix to diagnal position 
    for idx in range(min(di, dj)):

        T = new_X[idx: , idx: ]
        i, j = np.unravel_index(T.argmax(), T.shape) # get the coords of the max element of T
        
        if threshold and T[i, j] < threshold:
            dm = idx # new_X[:dm, :dm] is done (0, 1, ..., dm-1) excluding dm
            break
        else:
            dm = idx+1 # new_X[:dm, :dm] will be done

        # swap row idx, idx+i
        tmp = new_X[idx, :].copy()
        new_X[idx, :] = new_X[idx+i, :].copy() 
        new_X[idx+i, :] = tmp 
        
        tmp = new_rows[idx]
        new_rows[idx] = new_rows[idx+i]
        new_rows[idx+i] = tmp

        # swap col idx, idx+j
        tmp = new_X[:, idx].copy()
        new_X[:, idx] = new_X[:, idx+j].copy() 
        new_X[:, idx+j] = tmp 
        
        tmp = new_cols[idx]
        new_cols[idx] = new_cols[idx+j]
        new_cols[idx+j] = tmp
        
    # 
    if dm == dj:
        pass
    elif dm < dj: # free columns

        col_dict = {}
        sorted_col_idx = np.arange(dm)
        free_col_idx = np.arange(dm, dj)
        linked_rowcol_idx = new_X[:, dm:].argmax(axis=0)
        
        for col in sorted_col_idx:
            col_dict[col] = [col]
        for col, key in zip(free_col_idx, linked_rowcol_idx): 
            if key in col_dict.keys():
                col_dict[key] = col_dict[key] + [col]
            else:
                col_dict[key] = [col]
                
            
        new_col_order = np.hstack([col_dict[key] for key in sorted(col_dict.keys())])
        
        # update new_X new_cols
        new_X = new_X[:, new_col_order].copy()
        new_cols = new_cols[new_col_order]
    else:
        raise ValueError("Unexpected situation: dm > dj")
    
    if transposed:
        new_X = new_X.T
        new_rows, new_cols = new_cols, new_rows
    return new_X, new_rows, new_cols 

def diag_matrix_rows(X):
    """Diagonalize a matrix as much as possible by only rearrange rows
    """
    di, dj = X.shape
    
    new_X = np.array(X.copy())
    new_rows = np.arange(di) 
    new_cols = np.arange(dj) 
    
    # free to move rows
    row_dict = {}
    free_row_idx = np.arange(di)
    linked_rowcol_idx = new_X.argmax(axis=1) # the column with max value for each row
    
    for row, key in zip(free_row_idx, linked_rowcol_idx): 
        if key in row_dict.keys():
            row_dict[key] = row_dict[key] + [row]
        else:
            row_dict[key] = [row]
            
    new_row_order = np.hstack([row_dict[key] for key in sorted(row_dict.keys())])
    # update new_X new_cols
    new_X = new_X[new_row_order, :].copy()
    new_rows = new_rows[new_row_order]
    
    return new_X, new_rows, new_cols 

def csv_to_h5ad(
    csv_file, h5ad_file, sep,
    overwrite=False,
    ):
    """
    """
    import anndata

    if not overwrite and os.path.isfile(h5ad_file):
        logging.info("skipped {}, already exists".format(h5ad_file))
        return 

    # read
    logging.info("reading in {}".format(csv_file))
    df = pd.read_csv(csv_file, sep=sep, index_col=0)

    # transform
    h5ad_mat = anndata.AnnData(df.values, 
                              obs=pd.DataFrame(index=df.index.values),
                              var=pd.DataFrame(index=df.columns.values),
                             )

    # write
    h5ad_mat.write(h5ad_file, compression='gzip')
    logging.info("saved {}".format(h5ad_file))
    return 

def merge_h5ad(mats):
    merged = mats[0].concatenate(*mats[1:])
    return merged

def h5ad2pd(mat_h5ad):
    return pd.DataFrame(mat_h5ad.X, index=mat_h5ad.obs.index, columns=mat_h5ad.var.index)

def get_mcc(df_mc, df_c, 
            base_call_cutoff=100, 
            sufficient_coverage_fraction=1, 
            drop_features=True,
            fillna=True, 
            check=True):
    """
    Arguments:
        - feature by sample matrices

    Output:
        - mcc matrix (feature by sample)
        - return np.nan for entries less than `base_call_cutoff`
            drop_features:
            - further remove features has less than `sufficient_coverage_fraction` fraction of samples passing `base_call_cutoff`
                fillna:
                - further fillna for nan entries by the mean of non-nan entries across samples
    """
    # the two dataframes must have the same dimensions
    if check:
        assert df_c.shape == df_mc.shape
        assert np.all(df_c.index.values == df_mc.index.values)
        assert np.all(df_c.columns.values == df_mc.columns.values)
        assert np.all((df_c.values-df_mc.values) >= 0)
    
    _, n = df_c.shape
    
    # get mcc matrix with kept bins and nan values for low coverage sites
    df_c_nan = df_c.copy()
    df_c_nan[df_c < base_call_cutoff] = np.nan
    df_mcc = df_mc/df_c_nan 
    
    if drop_features:
        # a feature (region/gene) is sufficiently covered in % of cells 
        condition = (df_c > base_call_cutoff).sum(axis=1) >= sufficient_coverage_fraction*n
        logging.info("Matrix size before pruning... "+ str(df_mcc.shape))
        df_mcc = df_mcc.loc[condition]
        logging.info("Matrix size after pruning... "+ str(df_mcc.shape))

        # imputation (missing value -> mean value of all cells)
        if fillna:
            logging.info('Imputing data... (No effect if sufficient_coverage_fraction=1)')
            means = df_mcc.mean(axis=1)
            fill_value = pd.DataFrame({col: means for col in df_mcc.columns})
            df_mcc.fillna(fill_value, inplace=True)
    
    return df_mcc

def dedup_array_elements(x, empty_string=''):
    """Replacing repeats with empty_string
    """
    newx = np.empty_like(x)
    newx[0] = x[0]
    for i in range(1, len(x)):
        if x[i-1] == x[i]:
            newx[i] = empty_string
        else:
            newx[i] = x[i]
    return newx

def logcpm(counts):
    """
    Args:
        - gene-cell matrix (pandas DataFrame)
    """
    cov = counts.sum(axis=0)
    logcpm = np.log10(counts.divide(cov, axis=1)*1000000 + 1)
    return logcpm

def logtpm(counts, gene_lengths):
    """
    Args:
        - gene-cell matrix (pandas DataFrame)
        - gene_lengths: a series indexed by gene_id
    """
    tpm = counts.divide(gene_lengths.loc[counts.index], axis=0)
    cov = tpm.sum(axis=0)
    logtpm = np.log10((tpm.divide(cov, axis=1))*1000000 + 1)
    return logtpm

def sparse_logcpm(gc_matrix, mode='logcpm', lib_size=[]):
    """
    """
    lib_size = np.array(lib_size)
    if np.size(lib_size) == 0:
        lib_size = gc_matrix.data.sum(axis=0)

    lib_size_inv = sparse.diags(np.ravel(1.0/(1e-7+lib_size)))
    cpm = (gc_matrix.data).dot(lib_size_inv*1e6).tocoo()

    if mode == 'logcpm':
        cpm.data = np.log10(cpm.data + 1)
    elif mode == 'cpm':
        pass

    gc_cpm = GC_matrix(
        gc_matrix.gene, 
        gc_matrix.cell, 
        cpm,
    )
    
    return gc_cpm

def sparse_logtpm(gc_matrix, gene_lengths):
    """
    gene_lengths: array like 
    
    """
    gene_lengths = np.array(gene_lengths)
    gene_length_inv = sparse.diags(np.ravel(1.0/gene_lengths))
    tmp = (gene_length_inv).dot(gc_matrix.data).tocoo()
    lib_size_inv = sparse.diags(np.ravel(1.0/tmp.sum(axis=0)))
    
    logtpm = tmp.dot(lib_size_inv*1e6).tocoo()
    logtpm.data = np.log10(logtpm.data + 1)

    gc_logtpm = GC_matrix(
        gc_matrix.gene, 
        gc_matrix.cell, 
        logtpm,
    )
    
    return gc_logtpm