import numpy as np
import pandas as pd
import os
import logging
from pebble import ProcessPool
import subprocess as sp

import utils

# merge files
# takes too long and no memory (187k cells ~7hrs, 20G after gzip)
def get_mC_matrices(
    cells, files, ftypes, output_format, 
    name_cols=['name'],
    names=[],
    strict=False,
    overwrite=False,
    CAREFUL=False, 
    sep='\t', 
    ):
    """each cell has a file (tsv table) with name, ftypes
    generate output files for each ftype
    output_format.format(ftype)
    """
    logging.info("processing {}".format(output_format))
    assert len(cells) == len(files)

    # if all of the ftypes for this batch exists, skip...
    if not overwrite:
        flag = 0
        for ftype in ftypes:
            file = output_format.format(ftype)
            if (
                os.path.isfile(file) 
                or os.path.isfile(file+'.gz') 
                or os.path.isfile(file+'.bgz')
               ):
                # exists
                flag += 1
        if flag == len(ftypes):
            logging.info("{} exists, skipped".format(output_format.format('*')))
            return 

    # open file handles
    fhs = {}
    for ftype in ftypes:
        file = output_format.format(ftype)
        fh = open(file, 'w')
        fhs[ftype] = fh

    # go over cell by cell 
    names = np.array(names)
    i = 0
    for cell_id, file in zip(cells, files):
        try:
            tbl = pd.read_csv(file, sep=sep,)
        except:
            if strict:
                raise ValueError("Unable to open "+file)
            else:
                logging.info("skip "+file+'...')
                continue
        i += 1 # the i th number of file readed in  

        tmp_names = tbl[name_cols].apply(lambda x: "_".join(x.astype(str)), axis=1).values
        if i == 1:
            # establish index
            if len(names) == 0:
                names = tmp_names
            # write header
            for ftype in ftypes:
                towrite = "cell_id\t"+"\t".join(names.tolist())+"\n"
                fhs[ftype].write(towrite)
        if CAREFUL and (~np.all(names == tmp_names)):
            # sort table according to names (first file or specified)
            tbl['_name'] = tmp_names
            tbl = tbl.set_index('_name').reindex(names, fill_value=0).reset_index().drop('_name', axis=1)

        # write one row in each file
        for ftype in ftypes:
            towrite = cell_id+"\t"+"\t".join(tbl[ftype].astype(str).tolist())+"\n"
            fhs[ftype].write(towrite)

    # close file handles
    for ftype in ftypes:
        fhs[ftype].close()

    # gzip
    for ftype in ftypes:
        output = output_format.format(ftype)
        sp.run(['gzip', '-f', output])

    logging.info("{} done".format(output_format))
    return 

def runbatch_get_mC_matrices(
    cells, 
    files,
    datasets, 
    ftypes, 
    output_template, 
    name_cols=['name'],
    names=[],
    sep='\t', 
    overwrite=False, 
    CAREFUL=True,
    nprocs=1,
    ):
    """
    """
    assert len(cells) == len(files)
    assert len(cells) == len(datasets)

    logging.info(
        '''
        Number of cells: {}
        Number of datasets: {}
        Output template: {}
        ftypes: {}
        name_cols: {}
        names: {}
        overwrite: {}
        CAREFUL: {}
        nprocs: {}
        Number of expected files (#datasets * #ftypes): {}
        '''.format(
            len(cells), 
            len(np.unique(datasets)),
            output_template,
            ftypes,
            name_cols,
            len(names),
            overwrite,
            CAREFUL,
            nprocs,
            len(np.unique(datasets))*len(ftypes),
            )
    )

    metatable = pd.DataFrame()
    metatable['cell'] = cells
    metatable['dataset'] = datasets
    metatable['file'] = files
    with ProcessPool(max_workers=nprocs, max_tasks=10) as pool:
        for dataset, subtable in metatable.groupby('dataset'):
            output_format = output_template.format(dataset)
            subcells = subtable['cell'].values
            subfiles = subtable['file'].values
            future = pool.schedule(
                get_mC_matrices,
                args=(subcells, subfiles, ftypes, output_format), 
                kwargs=dict(
                    name_cols=name_cols, 
                    names=names,
                    overwrite=overwrite, 
                    CAREFUL=CAREFUL,
                    sep=sep,
                    ),
                ) 
            future.add_done_callback(utils.task_done)



# # example usage
# output_template = '../data/summary/bins100kb_{}_{{}}.tsv' # .format(dataset).format(context/ftype) 

# # get cells table
# cells_table = pd.read_csv('/cndd2/fangming/projects/m3c_human_dev/data/cells.tsv', sep='\t')
# cells_table['path_binc'] = cells_table['path_allc'].apply(
#     lambda x: os.path.join('../data/bins', 'bins100kb_'+os.path.basename(x))
# )
# print(cells_table.shape)
# cells_table.head()

# cells = cells_table['cell_id'].values
# datasets = cells_table['dataset'].values
# files = cells_table['path_binc'].values

# # run
# ftypes = ['CH', 'mCH', 'CG', 'mCG']
# name_cols = ['chr', 'bin']
# overwrite = False
# CAREFUL = False
# nprocs = 8

# # ### test !!!!
# cells_table_test = cells_table[cells_table['dataset']=='20210111_UMB1863'].head(10)
# cells = cells_table_test['cell_id']
# datasets = cells_table_test['dataset']
# files = cells_table_test['path_binc']
# name_cols = ['chr', 'bin']
# overwrite = True
# nprocs = 1

# batchrun_get_mC_matrices(
#     cells, 
#     files,
#     datasets, 
#     ftypes, 
#     output_template, 
#     name_cols=name_cols,
#     overwrite=overwrite, 
#     CAREFUL=CAREFUL,
#     nprocs=nprocs,
#     )
