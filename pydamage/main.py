#!/usr/bin/env python3
import pysam
from . import utils
import multiprocessing
from functools import partial
from statsmodels.stats.multitest import multipletests
from . import damage
import pandas as pd


def main(bam, wlen, show_al, mini, process, output, verbose):
    """Main function

    Args:
        bam (str): Path to alignment (sam/bam/cram) file
        wlen (int): window length
        show_al(bool): print alignments representations
        verbose(bool): verbose mode
    """
    mode = utils.check_extension(bam)
    alf = pysam.AlignmentFile(bam, mode)
    refs = list(alf.references)

    if len(refs) == 0:
        print(f"No aligned sequences in {bam}")
        return([])

    proc = min(len(refs), process)

    # Simple loop for debugging
    # all_res = []
    # for ref in refs:
    #     res = damage.test_ct(bam=bam, ref=ref, wlen=wlen,
    #                          min_al=mini, show_al=show_al,
    #                          mode=mode, process=process, verbose=verbose)
    #     if res:
    #         all_res.append(res)

    test_ct_partial = partial(damage.test_ct, bam=bam, wlen=wlen,
                              min_al=mini, show_al=show_al,
                              mode=mode, process=process, verbose=verbose)
    with multiprocessing.Pool(proc) as p:
        res = p.map(test_ct_partial, refs)
        filt_res = [i for i in res if i]
    df = pd.DataFrame(filt_res)
    multipletests(df['pvalue'], method='fdr_bh')
    df['qvalue'] = multipletests(df['pvalue'], method='fdr_bh')[1]
    df = df[['unif_pmin', 'geom_p', 'geom_pmin',
             'geom_pmax', 'pvalue', 'qvalue', 'reference']]
    df.sort_values(by=['qvalue'], inplace=True)
    df.to_csv(output)
