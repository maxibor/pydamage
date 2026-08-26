#!/usr/bin/env python3

import logging
import os
import shutil
import sys
from typing import Tuple

import numpy as np
import pandas as pd
from pysam import AlignmentFile
from statsmodels.stats.multitest import multipletests


def check_extension(filename: str) -> str:
    """Check alignment file format to give correct opening mode.

    Args:
        filename (str): Path to alignment file

    Returns:
        str: Opening mode ('rb', 'r', or 'rc')

    Raises:
        ValueError: If file extension is not supported
    """
    extension = filename.split(".")[-1].lower()
    modes = {"bam": "rb", "sam": "r", "cram": "rc"}
    if extension in modes:
        return modes[extension]
    raise ValueError(f"{extension} file extension not supported")


def makedir(dirpath: str, confirm: bool = True, force: bool = False):
    """Create directory, optionally confirming before overwrite."""
    if not dirpath:
        return
    if os.path.exists(dirpath):
        if confirm and not force:
            logging.warning(
                f"Result directory, {dirpath}, already exists, it will be overwritten"
            )
            if input("Do You Want To Continue? (y|n) ").lower() != "y":
                sys.exit(0)
        shutil.rmtree(dirpath)

    os.makedirs(dirpath, exist_ok=True)


def pandas_processing(res_dict: list, wlen: int, g2a: bool) -> pd.DataFrame:
    """Performs Pandas formatting and multiple testing correction on PyDamage results.

    Args:
        res_dict (list): List of result dictionaries from test_damage
        wlen (int): Window size
        g2a (bool): Output G to A transitions

    Returns:
        pd.DataFrame: Formatted DataFrame
    """
    df = pd.DataFrame(res_dict)
    if len(res_dict) == 0:
        return df

    valid_pvalues = df["pvalue"].dropna()
    if len(valid_pvalues) > 0:
        qvalues = pd.Series(
            multipletests(valid_pvalues, method="fdr_bh")[1],
            index=valid_pvalues.index,
            name="qvalue",
        )
        df = df.merge(qvalues, left_index=True, right_index=True, how="outer")
    else:
        df["qvalue"] = np.nan

    columns = [
        "p0",
        "p0_stdev",
        "p",
        "p_stdev",
        "pmin",
        "pmin_stdev",
        "pmax",
        "pmax_stdev",
        "pvalue",
        "qvalue",
        "RMSE",
        "reference",
        "nb_reads_aligned",
        "coverage",
        "reflen",
    ] + [f"CtoT-{i}" for i in range(wlen)]

    if g2a:
        columns += [f"GtoA-{i}" for i in range(wlen)]

    available_cols = [c for c in columns if c in df.columns]
    df = df[available_cols]
    df.rename(
        columns={
            "p0": "null_model_p0",
            "p0_stdev": "null_model_p0_stdev",
            "p": "damage_model_p",
            "p_stdev": "damage_model_p_stdev",
            "pmin": "damage_model_pmin",
            "pmin_stdev": "damage_model_pmin_stdev",
            "pmax": "damage_model_pmax",
            "pmax_stdev": "damage_model_pmax_stdev",
        },
        inplace=True,
    )
    if "qvalue" in df.columns:
        df.sort_values(by=["qvalue"], inplace=True)
    if "reference" in df.columns:
        df.set_index("reference", inplace=True)
    df.dropna(axis=1, how="all", inplace=True)
    df = df.round(3)
    return df


def df_to_csv(df: pd.DataFrame, outdir: str, outfile: str = "pydamage_results.csv"):
    """Write PyDamage results to disk."""
    df = df.round(3)
    if not outdir:
        outdir = "."
    outpath = os.path.join(outdir, outfile)
    df.to_csv(outpath)


def RMSE(residuals: np.ndarray) -> float:
    """Compute Root Mean Square Error."""
    return float(np.sqrt(np.mean(residuals**2)))


def prepare_bam(bam: str, minlen: int) -> Tuple[Tuple, str]:
    """Checks for file extension and returns tuple of mapped references >= minlen."""
    mode = check_extension(bam)
    alf = AlignmentFile(bam, mode)

    if not alf.has_index():
        logging.error(
            f"BAM file {bam} has no index. Sort BAM file and provide index "
            "before running pydamage."
        )
        alf.close()
        sys.exit(1)

    present_refs = set()
    for ref_stat, ref_len in zip(alf.get_index_statistics(), alf.lengths):
        refname = ref_stat[0]
        nb_mapped_reads = ref_stat[1]
        if nb_mapped_reads > 0 and ref_len >= minlen:
            present_refs.add(refname)
    alf.close()
    return tuple(sorted(present_refs)), mode
