#!/usr/bin/env python3
import pysam
from pydamage import utils
import multiprocessing
from functools import partial
from pydamage import damage
from pydamage.plot import damageplot
from pydamage.exceptions import PyDamageWarning
from pydamage.accuracy_model import prepare_data, load_model, fit_model
import sys
from tqdm import tqdm
import warnings
from pydamage import __version__


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
    ##########################
    ##########################

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
        filt_res = [
            damage.test_damage(
                ref=None,
                bam=bam,
                mode=mode,
                wlen=wlen,
                show_al=show_al,
                process=process,
                verbose=verbose,
            )
        ]
    else:
        with multiprocessing.Pool(proc) as p:
            res = list(tqdm(p.imap(test_damage_partial, refs), total=len(refs)))
        filt_res = [i for i in res if i]

    print(f"{len(filt_res)} contig(s) successfully analyzed by Pydamage")
    if len(filt_res) == 0:
        warnings.warn(
            "No alignments were found, check your alignment file", PyDamageWarning
        )

    if plot and len(filt_res) > 0:
        print("\nGenerating Pydamage plots")
        plotdir = f"{outdir}/plots"
        utils.makedir(plotdir, confirm=False)

        plot_partial = partial(damageplot, outdir=plotdir, wlen=wlen)
        with multiprocessing.Pool(proc) as p:
            list(tqdm(p.imap(plot_partial, filt_res), total=len(filt_res)))
    df_pydamage = utils.pandas_processing(res_dict=filt_res, wlen=wlen)

    acc_model = load_model()
    prep_df_glm = prepare_data(df_pydamage)
    df_glm = fit_model(prep_df_glm, acc_model)

    df = df_glm.merge(df_pydamage, left_index=True, right_index=True)

    utils.df_to_csv(df, outdir)
    return df
