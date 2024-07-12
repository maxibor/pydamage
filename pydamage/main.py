#!/usr/bin/env python3
import pysam
from pydamage import utils
import multiprocessing
from functools import partial
from pydamage import damage
from pydamage.plot import damageplot
from pydamage.exceptions import PyDamageWarning, AlignmentFileError
from pydamage.accuracy_model import prepare_data, glm_predict
from pydamage.rescale import rescale_bam
from pydamage.models import glm_model_params
import os
import sys
from tqdm import tqdm
import warnings
from pydamage import __version__
from collections import ChainMap

# def pydamage_analyze(
#     bam,
#     wlen=30,
#     show_al=False,
#     process=1,
#     outdir="",
#     plot=False,
#     verbose=False,
#     force=False,
#     group=False,
# ):

#     if group:
#         pydamage_analyze_group(
#             bam, wlen, show_al, process, outdir, plot, verbose, force
#         )

#     else:
#         pydamage_analyze_multi(
#             bam, wlen, show_al, process, outdir, plot, verbose, force
#         )


def pydamage_analyze(
    bam,
    wlen=30,
    minlen=0,
    show_al=False,
    process=1,
    outdir="",
    plot=False,
    verbose=False,
    force=False,
    group=False,
):
    """Runs the pydamage analysis for each reference separately

    Args:
        bam(str): Path to alignment (sam/bam/cram) file
        minlen(int): minimum reference sequence length threshold
        wlen(int): window length
        show_al(bool): print alignments representations
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

    refs, mode = utils.prepare_bam(bam, minlen=minlen)

    proc = max(1, min(len(refs), process))

    ##########################
    # Simple loop for debugging
    ##########################
    # filt_res = []
    # for ref in refs:
    #     res = damage.test_damage(
    #         bam=bam,
    #         ref=ref,
    #         wlen=wlen,
    #         show_al=show_al,
    #         mode=mode,
    #         process=process,
    #         verbose=verbose,
    #     )
    #     if res:
    #         filt_res.append(res)
    #     break
    ######################
    # Multiprocessing code
    ######################

    test_damage_partial = partial(
        damage.test_damage,
        bam=bam,
        mode=mode,
        wlen=wlen,
        show_al=show_al,
        process=process,
        verbose=verbose,
    )
    print("Estimating and testing Damage")
    if group:
        filt_res, read_dict = damage.test_damage(
            ref=None,
            bam=bam,
            mode=mode,
            wlen=wlen,
            show_al=show_al,
            process=process,
            verbose=verbose,
        )
        filt_res = [filt_res]
    else:
        if len(refs) > 0:
            with multiprocessing.Pool(proc) as p:
                res = list(tqdm(p.imap(test_damage_partial, refs), total=len(refs)))
            filt_res, read_dicts = zip(*res)
            filt_res = [i for i in filt_res if i]
            read_dicts = [i for i in read_dicts if i]
            read_dict = dict(ChainMap(*read_dicts))
        else:
            raise AlignmentFileError(
                "No reference sequences were found in alignment file"
            )

    ######################
    ######################

    print(f"{len(filt_res)} contig(s) analyzed by Pydamage")
    if len(filt_res) == 0:
        warnings.warn(
            "No alignments were found, check your alignment file", PyDamageWarning
        )

    rescale = True
    if rescale:
        print("\nRescaling quality scores")
        rescale_bam(
            bam=bam,
            threshold=0.5,
            alpha=0.05,
            damage_dict=filt_res,
            read_dict=read_dict,
            outname=os.path.join(outdir, "rescaled.bam"),
        )

    if plot and len(filt_res) > 0:
        print("\nGenerating Pydamage plots")
        plotdir = f"{outdir}/plots"
        utils.makedir(plotdir, confirm=False)

        plot_partial = partial(damageplot, outdir=plotdir, wlen=wlen)
        with multiprocessing.Pool(proc) as p:
            list(tqdm(p.imap(plot_partial, filt_res), total=len(filt_res)))
    df_pydamage = utils.pandas_processing(res_dict=filt_res, wlen=wlen)

    prep_df_glm = prepare_data(df_pydamage)
    df_glm = glm_predict(prep_df_glm, glm_model_params)

    df = df_glm.merge(df_pydamage, left_index=True, right_index=True)

    utils.df_to_csv(df, outdir)
    return df
