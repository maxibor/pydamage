#!/usr/bin/env python3
import pysam
from pydamage import utils
import multiprocessing
from functools import partial
from pydamage import damage
from pydamage.plot import damageplot
import pandas as pd
import sys
from tqdm import tqdm
from pydamage import __version__

def analyze(bam, wlen=30, show_al=False, mini=2000, cov=0.5, process=1, outdir="", plot = False, verbose=False, force=False):
    """Runs the pydamage analysis

    Args:
        bam(str): Path to alignment (sam/bam/cram) file
        wlen(int): window length
        show_al(bool): print alignments representations
        mini(int):  Minimum numbers of reads aligned  to consider contigs
        cov(float): Minimum coverage to consider contig
        process(int):  Number of  processes for parellel computing
        outdir(str): Path to output directory
        verbose(bool): verbose mode
        force(bool): force overwriting of results directory
    Returns:
        pd.DataFrame: pandas DataFrame containg pydamage results

    """
    if verbose:
        print(f"Pydamage version {__version__}\n")
    utils.makedir(outdir, force=force)

    mode = utils.check_extension(bam)
    alf = pysam.AlignmentFile(bam, mode)

    if not alf.has_index():
        print(f"BAM file {bam} has no index. Sort BAM file and provide index "
               "before running pydamage.")
        sys.exit(1)

    refs = list(alf.references)

    if len(refs) == 0:
        print(f"No aligned sequences in {bam}")
        return([])

    proc = min(len(refs), process)

    # Simple loop for debugging
    # all_res = []
    # for ref in refs:
    #     res = damage.test_damage(bam=bam, ref=ref, wlen=wlen,
    #                          min_al=mini, show_al=show_al,
    #                          mode=mode, process=process, verbose=verbose)
    #     if res:
    #         all_res.append(res)
    
    # print(all_res)

    test_damage_partial = partial(damage.test_damage, bam=bam, wlen=wlen,
                              min_al=mini, min_cov=cov, show_al=show_al,
                              mode=mode, process=process, verbose=verbose)
    with multiprocessing.Pool(proc) as p:
        res = p.map(test_damage_partial, refs)
        filt_res = [i for i in res if i]
    if plot and len(filt_res) > 0:
        print("\nGenerating pydamage plots")
        for ref in tqdm(filt_res):
            dam_plot = damageplot(damage_dict=ref, wlen=wlen, qlen = ref['qlen'], outdir=outdir)
            dam_plot.makedir()
            dam_plot.draw()
    df = utils.pandas_processing(res_dict=filt_res, outdir=outdir)
    return(df)
    
