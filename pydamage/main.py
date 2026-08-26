#!/usr/bin/env python3

from collections import ChainMap
from functools import partial
import logging
import multiprocessing
import os
import warnings
import pandas as pd
import pysam
from tqdm import tqdm

from pydamage import __version__
from pydamage import utils
from pydamage import damage
from pydamage.exceptions import PyDamageWarning, AlignmentFileError
from pydamage.models import prepare_data, glm_predict, glm_model_params
from pydamage.rescale import rescale_bam


def get_mp_context():
    """Get the most efficient multiprocessing context available on the platform."""
    available = multiprocessing.get_all_start_methods()
    if "fork" in available:
        return multiprocessing.get_context("fork")
    return multiprocessing.get_context()


def pydamage_analyze(
    bam: str,
    wlen: int = 30,
    minlen: int = 0,
    threshold: float = 0.5,
    show_al: bool = False,
    process: int = 1,
    outdir: str = "",
    plot: bool = False,
    verbose: bool = False,
    force: bool = False,
    group: bool = False,
    rescale: bool = False,
    subsample: float = None,
    no_ga: bool = False,
) -> pd.DataFrame:
    """Runs the PyDamage analysis for each reference contig or whole file.

    Args:
        bam (str): Path to alignment (SAM/BAM/CRAM) file
        wlen (int): Window length for damage estimation
        minlen (int): Minimum reference sequence length threshold
        threshold (float): Predicted accuracy threshold
        show_al (bool): Print alignment representations
        process (int): Number of processes for parallel computing
        outdir (str): Path to output directory
        plot (bool): Write damage fitting plots to disk
        verbose (bool): Verbose logging mode
        force (bool): Force overwriting of results directory
        group (bool): Treat entire BAM file as single reference
        rescale (bool): Rescale base quality scores using the PyDamage damage model
        subsample (float): Subsample fraction of reads (0-1)
        no_ga (bool): Do not use G->A transitions

    Returns:
        pd.DataFrame: DataFrame containing PyDamage results
    """
    g2a = not no_ga
    if subsample and rescale:
        raise ValueError("Cannot use subsample and rescale together")

    if no_ga and rescale:
        logging.warning("G to A transitions will not be used for rescaling!")

    if verbose:
        logging.info(f"Pydamage version {__version__}\n")

    utils.makedir(outdir, force=force)
    refs, mode = utils.prepare_bam(bam, minlen=minlen)
    proc = max(1, min(len(refs), process))

    if group:
        if verbose:
            logging.info("Estimating and testing Damage")
        filt_res_item, read_dict = damage.test_damage(
            ref=None,
            bam=bam,
            mode=mode,
            wlen=wlen,
            g2a=g2a,
            subsample=subsample,
            show_al=show_al,
            process=1,
            verbose=verbose,
        )
        filt_res = [filt_res_item] if filt_res_item else []
    else:
        if len(refs) > 0:
            if proc == 1:
                # Single-process fast-path: single open AlignmentFile handle
                al_handle = pysam.AlignmentFile(bam, mode=mode, threads=1)
                res = []
                for ref in tqdm(refs, desc="Estimating and testing Damage"):
                    r = damage.test_damage_with_handle(
                        ref=ref,
                        al_handle=al_handle,
                        wlen=wlen,
                        g2a=g2a,
                        subsample=subsample,
                        show_al=show_al,
                        verbose=verbose,
                        need_read_dict=rescale,
                    )
                    res.append(r)
                al_handle.close()
            else:
                # Multiprocess worker pool with persistent handles and optimal start context
                ctx = get_mp_context()
                chunksize = max(1, len(refs) // (proc * 4))
                worker_func = partial(
                    damage.test_damage_worker_task,
                    wlen=wlen,
                    g2a=g2a,
                    subsample=subsample,
                    show_al=show_al,
                    verbose=verbose,
                    need_read_dict=rescale,
                )
                with ctx.Pool(
                    proc,
                    initializer=damage.init_worker_handle,
                    initargs=(bam, mode),
                ) as p:
                    res = list(
                        tqdm(
                            p.imap(worker_func, refs, chunksize=chunksize),
                            total=len(refs),
                            desc="Estimating and testing Damage",
                        )
                    )

            filt_res_unfiltered, read_dicts_unfiltered = zip(*res)
            filt_res = [i for i in filt_res_unfiltered if i]
            if rescale:
                read_dicts = [i for i in read_dicts_unfiltered if i]
                read_dict = dict(ChainMap(*read_dicts)) if read_dicts else {}
            else:
                read_dict = {}
        else:
            raise AlignmentFileError("No reference sequences were found in alignment file")

    if not group and verbose:
        logging.info(f"{len(filt_res)} contig(s) analyzed by Pydamage")

    if len(filt_res) == 0:
        warnings.warn("No alignments were found, check your alignment file", PyDamageWarning)
        df_pydamage = pd.DataFrame()
        return df_pydamage

    if plot and len(filt_res) > 0:
        from pydamage.plot import damageplot
        if verbose:
            logging.info("Generating Pydamage plots")
        plotdir = os.path.join(outdir, "plots") if outdir else "plots"
        utils.makedir(plotdir, confirm=False)

        plot_partial = partial(damageplot, outdir=plotdir, wlen=wlen, plot_g2a=g2a)
        if proc == 1:
            for item in tqdm(filt_res, desc="Generating plots"):
                plot_partial(item)
        else:
            ctx = get_mp_context()
            chunksize = max(1, len(filt_res) // (proc * 4))
            with ctx.Pool(proc) as p:
                list(
                    tqdm(
                        p.imap(plot_partial, filt_res, chunksize=chunksize),
                        total=len(filt_res),
                        desc="Generating plots",
                    )
                )

    df_pydamage = utils.pandas_processing(res_dict=filt_res, wlen=wlen, g2a=g2a)
    prep_df_glm = prepare_data(df_pydamage)
    df_glm = glm_predict(prep_df_glm, glm_model_params)

    df = df_glm.merge(df_pydamage, left_index=True, right_index=True)
    utils.df_to_csv(df, outdir)

    if rescale:
        rescale_bam(
            bam=bam,
            threshold=threshold,
            alpha=0.05,
            wlen=wlen,
            damage_dict=df.to_dict(),
            read_dict=read_dict,
            grouped=group,
            outname=os.path.join(outdir, "pydamage_rescaled.bam") if outdir else "pydamage_rescaled.bam",
            threads=process,
        )

    return df
