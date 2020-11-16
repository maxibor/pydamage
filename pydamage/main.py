#!/usr/bin/env python3
import pysam
from pydamage import utils
import multiprocessing
from functools import partial
from pydamage import damage
from pydamage.plot import damageplot
from pydamage.exceptions import AlignmentFileError
from pydamage.accuracy_model import prepare_data, load_model, fit_model
import sys
from tqdm import tqdm
import warnings
from pydamage import __version__


def analyze(
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

    if group:
        analyze_group(bam,
                      wlen,
                      show_al,
                      process,
                      outdir,
                      plot,
                      verbose,
                      force)

    else:
        analyze_multi(bam,
                      wlen,
                      show_al,
                      process,
                      outdir,
                      plot,
                      verbose,
                      force)


def analyze_multi(
    bam,
    wlen=30,
    show_al=False,
    process=1,
    outdir="",
    plot=False,
    verbose=False,
    force=False,
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
        wlen=wlen,
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
    if len(filt_res) == 0 :
        raise AlignmentFileError("Check your alignment file")

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

    df = df_glm.merge(df_pydamage, left_index=True, right_index=True)

    utils.df_to_csv(df, outdir)
    return df


def analyze_group(
    bam,
    wlen=30,
    show_al=False,
    process=1,
    outdir="",
    plot=False,
    verbose=False,
    force=False,
):
    """Runs the pydamage analysis with all references grouped as one

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

    get_damage_group_partial = partial(
        damage.get_damage_group,
        bam=bam,
        wlen=wlen,
        show_al=show_al,
        mode=mode,
        process=process,
    )
    print("Estimating and testing Damage")
    with multiprocessing.Pool(proc) as p:
        res = list(
            tqdm(p.imap(get_damage_group_partial, refs), total=len(refs)))
    ct_data = []
    ga_data = []
    cc_data = []
    all_bases = []
    cov = 0
    nb_ref = 0
    nb_reads_aligned = 0
    reflen = 0
    for i in res:
        ct_data += i[0]
        ga_data += i[1]
        cc_data += i[2]
        all_bases += i[5]
        cov += i[6]
        nb_ref += 1
        nb_reads_aligned += i[7]
        reflen += i[8]
    cov = cov / nb_ref

    if nb_reads_aligned == 0:
        raise AlignmentFileError("No Alignments were found\nCheck your alignment file")


    damage_dict = damage.test_damage_group(
        ct_data,
        ga_data,
        cc_data,
        all_bases,
        nb_reads_aligned,
        cov,
        reflen,
        wlen,
        verbose,
    )

    if plot:
        print("\nGenerating Pydamage plot")
        plotdir = f"{outdir}/plots"
        utils.makedir(plotdir, confirm=False)
        damageplot(damage_dict, outdir=plotdir)

    df_pydamage = utils.pandas_group_processing(res_dict=damage_dict)

    acc_model = load_model()
    prep_df_glm = prepare_data(df_pydamage)
    df_glm = fit_model(prep_df_glm, acc_model)

    df = df_glm.merge(df_pydamage, left_index=True, right_index=True)

    utils.df_to_csv(df, outdir)
    return df
