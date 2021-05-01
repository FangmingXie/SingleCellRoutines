#!/usr/bin/env python3
DESCRIPTION = '''
A python wrapper of a bash/awk worker that counts mC level for each of the 
tri-nucleotide contexts
'''

from __init__scr import *
import argparse
import subprocess as sp
import numpy as np
import time
from pebble import ProcessPool
import logging
import os

import utils

def allc_count_context_worker_wrap(
    input_allc_file, 
    output_file,
    compress=True,
    overwrite=False,
    dirname=DIRNAME, # package directory
    ):
    """
    """
    logging.info("processing: {}".format(input_allc_file))
    if not overwrite:
        if os.path.isfile(output_file) or os.path.isfile(output_file+'.gz') or os.path.isfile(output_file+'.bgz'):
            logging.info("File exists "+output_file+", skipping...")
            return 0
    sp.run([os.path.join(dirname, "_allc_count_contexts_worker.sh"), input_allc_file, output_file])
    if compress:
        utils.compress(output_file)
    logging.info("Done. Results saved to {}".format(output_file))
    return 

def run_allc_count_contexts(
    input_allc_files, 
    output_prefix, 
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
                """.format(nprocs, len(input_allc_files)))

    output_files = [
        output_prefix+"_"+os.path.basename(input_allc_file).replace('.tsv.gz', '.tsv')
        for input_allc_file in input_allc_files] 

    output_dir = os.path.dirname(output_prefix) 
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # parallelized processing
    with ProcessPool(max_workers=nprocs, max_tasks=10) as pool:
        for input_allc_file, output_file in zip(input_allc_files, output_files):
            future = pool.schedule(allc_count_context_worker_wrap, 
                                   args=(input_allc_file, output_file,), 
                                   kwargs={
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
        "-o", "--output_prefix", 
        required=True,
        help="output directory and file prefix;\
              outputs will be named as output_prefix+'_'+input_name[.tsv.gz]",
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

if __name__  == '__main__':
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
        
    output_prefix = args.output_prefix

    compress = args.compress
    overwrite = args.overwrite
    nprocs = args.nprocs
    timeout = args.timeout

    logging.info(
        """ mC genomewide non-overlapping bins counting:
            Allc tables: {}
            Output prefix: {}
            Compress: {}
            Overwrite: {}
            Number of processes: {}
            Timeout: {}
        """.format(
            len(input_allc_files), 
            output_prefix, 
            compress,
            overwrite,
            nprocs, 
            timeout,
    ))
    
    run_allc_count_contexts(
        input_allc_files, 
        output_prefix, 
        compress=compress, 
        overwrite=overwrite,
        nprocs=nprocs,
        timeout=timeout,
    )
