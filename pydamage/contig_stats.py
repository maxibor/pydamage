#!/usr/bin/env python3

from io import StringIO
import numpy as np
import pandas as pd
import pysam


def compute_coverage_sum(pysam_cov) -> int:
    """Computes total base coverage sum of a reference from count_coverage output.

    Args:
        pysam_cov (tuple or np.ndarray): Four tuples/arrays of base coverage (A, C, G, T)
    Returns:
        int: Total sum of coverage across all bases
    """
    if isinstance(pysam_cov, np.ndarray):
        return int(np.sum(pysam_cov))
    return int(sum(sum(c) for c in pysam_cov))


def get_contig_stats(al_handle: pysam.AlignmentFile, ref: str) -> tuple:
    """Get contig length, average coverage, and aligned read count."""
    reflen = al_handle.get_reference_length(ref)
    cov_counts = al_handle.count_coverage(contig=ref)
    cov_sum = sum(sum(c) for c in cov_counts)
    nb_reads_aligned = al_handle.count(contig=ref)
    cov = cov_sum / reflen if reflen > 0 else 0.0
    return reflen, cov, nb_reads_aligned


def get_ref_stats(bam: str) -> tuple:
    """Get summary statistics for the entire BAM file."""
    cov_stats = pysam.coverage(bam)
    df = pd.read_table(StringIO(cov_stats))
    df = df.query("numreads > 0")
    total_endpos = df["endpos"].sum()
    if total_endpos > 0:
        cov = float((df["meandepth"] * df["endpos"]).sum() / total_endpos)
    else:
        cov = 0.0
    reflen = int(total_endpos)
    nb_reads_aligned = int(df["numreads"].sum())
    return reflen, cov, nb_reads_aligned