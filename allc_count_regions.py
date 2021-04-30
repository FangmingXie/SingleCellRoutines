#!/usr/bin/env python3
"""
"""
DESCRIPTION = '''
Read in allc tables (with tabix open) and generate summarized mc levels for given regions. 
Inputs: allctables, a bedfile (a list of regions), a list of contexts.
Outputs: one file for each allc table for a list of context.
'''
from __init__scr import *

import argparse
import tabix
from pebble import ProcessPool
import logging
import os
import pandas as pd

import utils

def mc_region_level_worker(
    allc_file, output_file, bed_file, 
    bed_file_name_column=False,
    contexts=CONTEXTS, 
    compress=True, 
    cap=2, 
    overwrite=False,
    ):
    """
    allc_file 
    bed file:
    """
    logging.info("Begin mc_region_level processing: {}".format(allc_file))
    # check overlap
    if not overwrite:
        if os.path.isfile(output_file) or os.path.isfile(output_file+'.gz') or os.path.isfile(output_file+'.bgz'):
            logging.info("File exists "+output_file+", skipping...")
            return 0

    # read in bed file and remove chr prefix
    if bed_file_name_column:
        columns = ['chr', 'start', 'end', 'name']
        df_gtf = pd.read_csv(bed_file, sep="\t", header=None, 
                names=columns, usecols=[0, 1, 2, 3], dtype={'chr': object, 'name': object})
    else:
        columns = ['chr', 'start', 'end']
        df_gtf = pd.read_csv(bed_file, sep="\t", header=None, 
                names=columns, usecols=[0, 1, 2], dtype={'chr': object})
    df_gtf.chr = df_gtf.chr.apply(lambda x: x[len('chr'):] if x.startswith('chr') else x)

    # auto detect whether the allc table has chr prefix or not
    allc_peek = utils.read_allc(allc_file, 
                                pindex=False, 
                                remove_chr=False,
                                nrows=1, 
                                )
    chrom_format = allc_peek.iloc[0]['chr']
    if chrom_format.startswith('chr'):
        allc_chr_prefix = True
    else:
        allc_chr_prefix = False

    # output header
    outfile = open(output_file, "w")
    for context in contexts:
        columns += ['m{}'.format(context), context]
    outfile.write('\t'.join(columns)+'\n')
    # write results region by region
    allc = tabix.open(allc_file)
    for i, row in df_gtf.iterrows():
        row_out = [str(row.chr), str(row.start), str(row.end)]
        if bed_file_name_column:
            row_out += [str(row['name'])]
        if not row_out[0].startswith('chr'): # enforce output has chr as prefix
            row_out[0] = 'chr'+row_out[0]

        if allc_chr_prefix: 
            records = list(allc.query('chr'+str(row['chr']), row['start'], row['end']))
        else:
            records = list(allc.query(row['chr'], row['start'], row['end']))

        for context in contexts:
            mc, c = utils.tabix_summary(records, context=context, cap=cap) # remove sites with total_c > 2
            row_out += [str(mc), str(c)] 

        outfile.write('\t'.join(row_out)+'\n')
    outfile.close()

    if compress:
        utils.compress(output_file, suffix='gz')
    logging.info("Done. Results saved to {}".format(output_file))
    return 0

def run_mc_region_level(
    input_allc_files, 
    input_bed_file, 
    output_prefix, 
    bed_file_name_column=False, 
    contexts=CONTEXTS,
    compress=True, 
    cap=2,
    overwrite=False,
    nprocs=1, 
    timeout=None,
    ):
    """
    run mc_gene_level in parallel
    """
    # assume certain structures in the inputs and outputs
    # allc_xxx.tsv.gz -> output_prefix + "_" + allc_xxx.tsv.gz
    # but the output_files should remove .gz suffix at first 
    output_files = [
        output_prefix+"_"+os.path.basename(input_allc_file).replace('.tsv.gz', '.tsv')
        for input_allc_file in input_allc_files] 

    nprocs = min(nprocs, len(input_allc_files))
    logging.info("""Begin run_mc_region_level.\n
                Number of processes:{}\n
                Number of allc_files:{}\n
                Bed file: {}\n
                """.format(nprocs, len(input_allc_files), input_bed_file))
    
    # parallelized processing
    with ProcessPool(max_workers=nprocs, max_tasks=10) as pool:
        for input_allc_file, output_file in zip(input_allc_files, output_files):
            future = pool.schedule(mc_region_level_worker, 
                                   args=(input_allc_file, output_file, input_bed_file), 
                                   kwargs={
                                        'bed_file_name_column': bed_file_name_column,
                                        'contexts': contexts, 
                                        'compress': compress,
                                        'cap': cap,
                                        'overwrite': overwrite,
                                        },
                                   timeout=timeout)
            future.add_done_callback(utils.task_done)
    # end parallel
    return

def create_parser():
    """

    """
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
        "-b", "--input_bed_file", 
        required=True,
        help="bed file",
    )
    parser.add_argument(
        "-bn", "--bed_file_name_column", 
        action='store_true',
        help="whether the bed file contains a name column",
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
        "--cap", 
        type=int, 
        default=2,
        help="Exclude totalC greater than this values;\
              turn it to 0 for bulk data",
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

    if isinstance(args.input_allc_files, list) and len(args.input_allc_files) > 0:
        input_allc_files = args.input_allc_files
    elif isinstance(args.input_allc_files_txt, str) and len(args.input_allc_files_txt) > 0:
        input_allc_files = utils.import_single_textcol(args.input_allc_files_txt)
    else:
        raise ValueError("no input files")
        
    input_bed_file = args.input_bed_file
    output_prefix = args.output_prefix

    bed_file_name_column = args.bed_file_name_column
    contexts = args.contexts
    compress = args.compress
    cap = args.cap
    overwrite = args.overwrite
    nprocs = args.nprocs
    timeout = args.timeout

    logging.info(
        """ CEMBA mc region level counting:
            Allc tables: {}
            Bed file: {} (Contain name col: {})
            Output prefix: {}
            Contexts: {}
            Compress: {}
            Cap: {}
            Overwrite: {}
            Number of processes: {}
            Timeout: {}
        """.format(
            len(input_allc_files), 
            input_bed_file, bed_file_name_column,
            output_prefix, 
            contexts,
            compress,
            cap, 
            overwrite,
            nprocs, 
            timeout,
    ))

    run_mc_region_level(
        input_allc_files, 
        input_bed_file,
        output_prefix, 
        bed_file_name_column=bed_file_name_column,
    	contexts=contexts,
        compress=compress,
        cap=cap,
    	overwrite=overwrite,
    	nprocs=nprocs,
        timeout=timeout,
        )