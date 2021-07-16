#!/usr/bin/env python
import glob
import logging
from pebble import ProcessPool
import utils

def runbatch_csv2h5ad(
    files,
    sep=',',
    overwrite=False,
    nprocs=1,
    ):
    """
    """
    logging.info(
        '''
        number of files: {}
        nprocs: {}
        '''.format(len(files), nprocs)
    )

    with ProcessPool(max_workers=nprocs, max_tasks=10) as pool:
        for file in files:
            output_file = file.replace('.tsv.gz', '.h5ad')
            future = pool.schedule(
                utils.csv_to_h5ad, 
                args=(file, output_file, sep),
                kwargs=dict(overwrite=overwrite),
                timeout=None,
                )
            future.add_done_callback(utils.task_done)
    return 
