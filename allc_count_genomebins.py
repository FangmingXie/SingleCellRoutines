#!/usr/bin/env python3
"""
"""
DESCRIPTION = '''
Read in allc tables and generate summarized mc levels over non-overlapping genomic bins. 
Inputs: allctables, a bin size, a genome size file, a list of contexts.
Outputs: one file for each allc table for a list of context.
'''

from __init__scr import *

import numpy as np
import pandas as pd
import time
import argparse
from collections import OrderedDict
from pebble import ProcessPool
import logging
import os
import pandas as pd
import subprocess as sp
import glob

import utils

def bin_allc_worker(
    allc_file, 
    genome_size_file,
    bin_size, 
    output_file,
    chromosomes=None, 
    contexts=CONTEXTS,
    compress=True,
    overwrite=False,
    ):

    logging.info("binc processing: {} {}".format(allc_file, contexts))

    if not overwrite:
        if os.path.isfile(output_file) or os.path.isfile(output_file+'.gz') or os.path.isfile(output_file+'.bgz'):
            logging.info("File exists "+output_file+", skipping...")
            return 0
    
    chrom_sizes = utils.get_chrom_lengths(genome_size_file)
    if chromosomes == None: 
        chromosomes = chrom_sizes.index.values

    # read allc
    df_allc = utils.read_allc(allc_file, pindex=False, compression='gzip')
    df_allc = df_allc.loc[df_allc.chr.isin(chromosomes)]

    # initiate
    chr_allc = np.array([])
    bins_allc = np.array([])
    mc_c_allc = OrderedDict()
    for context in contexts:
        mc_c_allc['m'+context] = np.array([])
        mc_c_allc[context] = np.array([])

    for chromosome, df in df_allc.groupby('chr'):
        # last bin (incomplete) is discarded
        bins = np.arange(0, chrom_sizes[chromosome], bin_size)

        # number of intervals is number of bin points -1 
        chrs = np.asarray([chromosome]*(len(bins)-1))
        bins_allc = np.concatenate([bins_allc, bins[:-1]])
        chr_allc = np.concatenate([chr_allc, chrs])

        for context in contexts:
            df_context = df.loc[df.context.isin(utils.get_expanded_context(context))]
            df_mc_c = df_context.groupby(pd.cut(df_context.pos, bins)).sum().fillna(0)[['mc', 'c']]
            mc_c_allc['m'+context] = np.concatenate([mc_c_allc['m'+context], df_mc_c.mc])
            mc_c_allc[context] = np.concatenate([mc_c_allc[context], df_mc_c.c])

    columns = ['chr', 'bin'] + [key for key in mc_c_allc]
    binc = pd.DataFrame(columns=columns)
    binc['chr'] = chr_allc.astype(object)
    binc['bin'] = bins_allc.astype(int)
    for key, value in mc_c_allc.items():
        binc[key] = value.astype(int) 

    binc.to_csv(output_file, na_rep='NA', sep="\t", header=True, index=False)
    # bgzip
    if compress:
        utils.compress(output_file, suffix='gz')
    logging.info("Done. Results saved to {}".format(output_file))

    return 

def run_bin_allc(
    input_allc_files, 
    genome_size_file,
    bin_size,
    output_prefix, 
    chromosomes=None,
    contexts=CONTEXTS, 
    compress=True, 
    overwrite=False,
    nprocs=1,
    timeout=None,
    ):
    """
    run bin_allc in parallel
    """
    # assume certain structures in the inputs and outputs
    # allc_xxx.tsv.gz -> output_prefix + "_" + allc_xxx.tsv.gz
    # but the output_files should remove .gz suffix at first 
    nprocs = min(nprocs, len(input_allc_files))
    logging.info("""Begin run bin allc.\n
                Number of processes:{}\n
                Number of allc_files:{}\n
                Genome size file: {}\n
                Bin size: {}\n
                """.format(nprocs, len(input_allc_files), genome_size_file, bin_size))

    output_files = [
        output_prefix+"_"+os.path.basename(input_allc_file).replace('.tsv.gz', '.tsv')
        for input_allc_file in input_allc_files] 

    output_dir = os.path.dirname(output_prefix) 
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # parallelized processing
    with ProcessPool(max_workers=nprocs, max_tasks=10) as pool:
        for input_allc_file, output_file in zip(input_allc_files, output_files):
            future = pool.schedule(bin_allc_worker, 
                                   args=(input_allc_file, genome_size_file, bin_size, output_file), 
                                   kwargs={
                                        'chromosomes': chromosomes,
                                        'contexts': contexts, 
                                        'compress': compress,
                                        'overwrite': overwrite,
                                        },
                                   timeout=timeout)
            future.add_done_callback(utils.task_done)
    # end parallel
    return 

def create_parser():
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    ) 

    parser.add_argument(
        "-i", "--input_allc_files", 
        nargs='+', 
        help="a list of paths to allc tables",
    )
    parser.add_argument(
        "-itxt", "--input_allc_files_txt", 
        type=str,
        help="a file containing a list of paths to allc tables",
    )
    parser.add_argument(
        "-g", "--genome_size_file", 
        type=str,
        required=True,
        help="chromosome sizes (two columns tsv)",
    )
    parser.add_argument(
        "-b", "--bin_size", 
        type=int,
        required=True,
        help="a bin size",
    )
    parser.add_argument(
        "-o", "--output_prefix", 
        required=True,
        help="output directory and file prefix;\
              outputs will be named as output_prefix+'_'+input_name[.tsv.gz]",
    )

    parser.add_argument(
        "-chr", "--chromosomes", 
        nargs='+', 
        default=None, 
        help="list of chromosomes",
    )
    parser.add_argument(
        "-c", "--contexts", 
        nargs='+', 
        default=CONTEXTS, 
        help="list of contexts: CH/CG/... default: CONTEXTS (CH CG CA)",
    )
    parser.add_argument(
        "-bgzip", "--compress", 
        action='store_true',
        help="bgzip the outputs",
    )
    parser.add_argument(
        "-n", "--nprocs", 
        type=int, 
        default=1,
        help="number of processes",
    )
    parser.add_argument(
        "-f", "--overwrite", 
        action='store_true',
        help="overwrite a file if it exists",
    )
    parser.add_argument(
        "-t", "--timeout", 
        type=int,
        default=None,
        help="? seconds per file time limit (timeout)",
    )
    return parser


if __name__ == '__main__':
    log = utils.create_logger()
    parser = create_parser()
    args = parser.parse_args()

    # input allc files
    if isinstance(args.input_allc_files, list) and len(args.input_allc_files) > 0:
        input_allc_files = args.input_allc_files
    elif isinstance(args.input_allc_files_txt, str) and len(args.input_allc_files_txt) > 0:
        input_allc_files = utils.import_single_textcol(args.input_allc_files_txt)
    else:
        raise ValueError("no input files")
        
    genome_size_file = args.genome_size_file
    bin_size = args.bin_size
    output_prefix = args.output_prefix

    chromosomes = args.chromosomes
    contexts = args.contexts
    compress = args.compress
    overwrite = args.overwrite
    nprocs = args.nprocs
    timeout = args.timeout

    logging.info(
        """ mC genomewide non-overlapping bins counting:
            Allc tables: {}
            Genome size file: {}
            Bin size: {}
            Output prefix: {}
            chromosomes: {}
            Contexts: {}
            Compress: {}
            Overwrite: {}
            Number of processes: {}
            Timeout: {}
        """.format(
            len(input_allc_files), 
            genome_size_file, 
            bin_size,
            output_prefix, 
            chromosomes,
            contexts,
            compress,
            overwrite,
            nprocs, 
            timeout,
    ))
    
    run_bin_allc(
        input_allc_files, 
        genome_size_file,
        bin_size,
        output_prefix, 
        chromosomes=chromosomes,
        contexts=contexts, 
        compress=compress, 
        overwrite=overwrite,
        nprocs=nprocs,
        timeout=timeout,
    )
