import numpy as np
import os
import sys
import shutil
from statsmodels.stats.multitest import multipletests
import pandas as pd



def check_extension(filename):
    extension = filename.split(".")[-1]
    modes = {'bam': 'rb', 'sam': 'r', 'cram': 'rc'}
    try:
        return(modes[extension])
    except KeyError:
        raise Exception(f"{extension} file extension not supported")


def get_x_y(data):
    xdata, counts = np.unique(data, return_counts=True)
    ydata = counts/counts.sum()
    return(xdata, ydata)


def makedir(dirpath, confirm=True, force=False):
    if os.path.exists(dirpath):
        if confirm and force==False:
            print(f"Result directory, {dirpath}, already exists, it will be overwritten")
            if input('Do You Want To Continue? (y|n) ').lower() != 'y':
                sys.exit()
        shutil.rmtree(dirpath)
    
    os.makedirs(dirpath)

def pandas_processing(res_dict, outdir):
    df = pd.DataFrame(res_dict)
    df['qvalue'] = multipletests(df['pvalue'], method='fdr_bh')[1]
    df = df[['unif_pmin', 'unif_pmin_stdev', 
             'geom_p', 'geom_p_stdev',
             'geom_pmin', 'geom_pmin_stdev',
             'geom_pmax', 'geom_pmax_stdev',
             'pvalue', 
             'qvalue', 
             'reference',
             'nb_reads_aligned',
             'coverage']+
             [f"CtoT-{i}" for i in range(df['qlen'].max())]+
             [f"GtoA-{i}" for i in range(df['qlen'].max())]]
    df.sort_values(by=['qvalue'], inplace=True)
    df.set_index("reference", inplace=True)
    df.dropna(axis=1, how='all', inplace=True)

    df.to_csv(f"{outdir}/pydamage_results.csv")
    return(df)