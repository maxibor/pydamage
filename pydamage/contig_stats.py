import numpy as np
from tqdm import tqdm
import pysam
from numba import jit
from functools import partial
import multiprocessing
import pandas as pd
from io import StringIO
import logging


@jit()
def compute_coverage_sum(pysam_cov):
    """Computes average coverage of a reference

    Args:
        pysam_cov (np.array): Four dimensional array of coverage for each base
        ref_len (int): length of reference
    Returns:
        float: mean coverage of reference
    """
    return (
        np.sum(
            np.sum(pysam_cov, axis=0)
        )
    )

def _get_contig_stats(ref, al_handle):
    """
    Get the contig stats for a given reference and contig.
    """
    reflen = al_handle.get_reference_length(ref)
    cov = compute_coverage_sum(
        np.array(al_handle.count_coverage(contig=ref))
    )
    nb_reads_aligned = al_handle.count(contig=ref)
    

    return reflen, cov, nb_reads_aligned

def get_contig_stats(al_handle, ref):
    reflen, cov_sum, nb_reads_aligned = _get_contig_stats(ref, al_handle)
    
    cov = cov_sum/reflen if reflen > 0 else 0
    
    return reflen, cov, nb_reads_aligned

def get_ref_stats(bam):
    cov_stats = pysam.coverage(bam)
    df = pd.read_table(StringIO(cov_stats))
    df = df.query("numreads > 0")
    cov = (df['meandepth']*df['endpos']).sum()/df['endpos'].sum()
    reflen = df['endpos'].sum()
    nb_reads_aligned = df['numreads'].sum()

    return reflen, cov, nb_reads_aligned



def get_ref_stats_pysam(al_handle):
    cov_sum = np.empty((0,1), dtype="uint32")
    nb_reads_aligned = np.empty((0,1), dtype="uint32")
    reflen = np.empty((0,1), dtype="uint32")

    for r in tqdm(al_handle.references):
        _reflen, _cov, _nb_reads_aligned = _get_contig_stats(r, al_handle)
        if _nb_reads_aligned > 0:
            nb_reads_aligned = np.append(nb_reads_aligned, _nb_reads_aligned)
            reflen = np.append(reflen, _reflen)
            cov_sum = np.append(cov_sum, _cov)
    
    cov = np.sum(cov_sum)/np.sum(reflen) if np.sum(reflen) > 0 else 0
    
    return np.sum(reflen), cov, np.sum(nb_reads_aligned)

def get_ref_stats_multi(bam, process):
    
    get_contig_stats_partial = partial(_get_contig_stats, bam=bam, process=process)
    al_handle = pysam.AlignmentFile(bam, mode="rb", threads=process)
    refs = al_handle.references
    al_handle.close()
    with multiprocessing.Pool(process) as p:
        res = np.array(p.map(get_contig_stats_partial, refs))
    return np.sum(res[:, 0]), np.sum(res[:, 1])/np.sum(res[:, 0]), np.sum(res[:, 2])