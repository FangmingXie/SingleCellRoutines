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