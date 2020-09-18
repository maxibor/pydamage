#!/usr/bin/env python3
import pysam
from pydamage import utils
import multiprocessing
from functools import partial
from pydamage import damage
from pydamage.plot import damageplot
from pydamage.accuracy_model import prepare_data, load_model, fit_model
import pandas as pd
import sys
from tqdm import tqdm
import warnings
from pydamage import __version__


def analyze(
    bam,
    fasta,
    wlen=30,
    show_al=False,
    mini=2000,
    cov=0.5,
    process=1,
    outdir="",
    plot=False,
    verbose=False,
    force=False,
):
    """Runs the pydamage analysis

    Args:
        bam(str): Path to alignment (sam/bam/cram) file
        fasta(str): Path to fasta file containing reference sequences
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

    if not verbose:
        warnings.filterwarnings("ignore")

    mode = utils.check_extension(bam)
    alf = pysam.AlignmentFile(bam, mode)

    if not alf.has_index():
        print(
            f"BAM file {bam} has no index. Sort BAM file and provide index "
            "before running pydamage."
        )
        sys.exit(1)

    refs = list(alf.references)

    if len(refs) == 0:
        print(f"No aligned sequences in {bam}")
        return []

    proc = min(len(refs), process)

    ##########################
    # Simple loop for debugging
    ##########################
    # filt_res = []
    # for ref in refs:
    #     res = damage.test_damage(bam=bam, ref=ref, wlen=wlen,
    #                              min_al=mini,  min_cov=cov, show_al=show_al,
    #                              mode=mode, process=process, verbose=verbose)
    #     if res:
    #         filt_res.append(res)
    ##########################
    ##########################

    test_damage_partial = partial(
        damage.test_damage,
        bam=bam,
        fasta=fasta,
        wlen=wlen,
        min_al=mini,
        min_cov=cov,
        show_al=show_al,
        mode=mode,
        process=process,
        verbose=verbose,
    )
    print("Estimating and testing Damage")
    with multiprocessing.Pool(proc) as p:
        res = list(tqdm(p.imap(test_damage_partial, refs), total=len(refs)))
    filt_res = [i for i in res if i]

    print(f"{len(filt_res)} contigs were successfully analyzed by Pydamage")

    if plot and len(filt_res) > 0:
        print("\nGenerating Pydamage plots")
        plotdir = f"{outdir}/plots"
        utils.makedir(plotdir, confirm=False)

        plot_partial = partial(damageplot, outdir=plotdir)
        with multiprocessing.Pool(proc) as p:
            list(tqdm(p.imap(plot_partial, filt_res), total=len(filt_res)))

    df_pydamage = utils.pandas_processing(res_dict=filt_res)

    acc_model = load_model()
    prep_df_glm = prepare_data(df_pydamage)
    df_glm = fit_model(prep_df_glm, acc_model)

    df = df_pydamage.merge(df_glm, left_index=True, right_index=True)

    utils.df_to_csv(df, outdir)
    return df
