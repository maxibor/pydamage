#!/usr/bin/env python3
import pysam
import utils
import multiprocessing
from functools import partial
import damage
import pandas as pd


def main(bam, wlen, show_al, mini, process, output):
    """Main function

    Args:
        bam (str): Path to alignment (sam/bam/cram) file
        wlen (int): window length
        show(bool): print alignments representations
    """
    mode = utils.check_extension(bam)
    alf = pysam.AlignmentFile(bam, mode)
    refs = list(alf.references)

    if len(refs) == 0:
        print(f"No aligned sequces in {bam}")
        return([])

    proc = min(len(refs), process)

    # Simple loop for debugging
    # all_res = []
    # for ref in refs:
    #     res = damage.test_ct(bam=bam, ref=ref, wlen=wlen,
    #                          min_al=mini, show_al=show_al,
    #                          mode=mode, process=process)
    #     if res:
    #         all_res.append(res)

    test_ct_partial = partial(damage.test_ct, bam=bam, wlen=wlen,
                              min_al=mini, show_al=show_al,
                              mode=mode, process=process)
    with multiprocessing.Pool(proc) as p:
        res = p.map(test_ct_partial, refs)
        filt_res = [i for i in res if i]
    print(pd.DataFrame(filt_res))
