import pysam
import numpy as np
from array import array
from pydamage.models import damage_model
from pydamage import __version__
from tqdm import tqdm
from numba import jit
import sys


@jit(
    cache=True,
    fastmath=True,
)
def phred_to_prob(qual):
    """Convert Phred quality score to probability

    Args:
        qual (array): Array of unsigned integer Phred quality scores

    Returns:
        np.array(int): Array of read error probabilities
    """
    return 10 ** (-qual / 10)


@jit(
    cache=True,
    fastmath=True,
)
def compute_new_prob(e, d):
    """Compute new probability of base calling  error accounting for ancient damage

    Args:
        e (np.array): Array of read error probabilities
        d (np.array): Array of damage probabilities
    Returns:
        np.array(int): Array of new read error probabilities
    """
    return np.round(-10 * np.log10(1 - np.multiply(1 - e, 1 - d)), 0).astype(np.int64)


def rescale_qual(read_qual, dmg_pmf, damage_bases, reverse):
    """Rescale quality scores using damage model
    Args:
        read_qual (array): Array of Phred quality scores
        dmg_pmf (array): Array of damage model probabilities
        damage_bases (array): Array of positions with damage
        reverse(bool): Read mapped to reverse strand
    Returns:
        np.array(int): Array of rescaled Phred quality scores
    """
    e = phred_to_prob(np.array(read_qual).astype(np.int64))
    if reverse:
        e = e[::-1]
    d = np.zeros(len(read_qual))
    d[damage_bases] = dmg_pmf[damage_bases]
    r = compute_new_prob(e, d)
    if reverse:
        r = r[::-1]
    return r


def rescale_bam(
    bam, threshold, alpha, wlen, damage_dict, read_dict, grouped, outname, threads
):
    """Rescale quality scores in BAM file using damage model

    Args:
        bam (str): Path to BAM file
        threshold (float): Predicted accuracy threshold
        alpha (float): Q-value threshold
        wlen (int): Window size for damage model
        damage_dict (dict): Damage model parameters
        read_dict (dict): Dictionary of read names
        grouped (bool): Grouped analysis
        outname (str): Path to output BAM file
        threads(int): Number of threads
    """
    with pysam.AlignmentFile(bam, "rb", threads=threads) as al:
        hd = al.header.to_dict()
        hd["PG"].append(
            {
                "ID": "pydamage",
                "PN": "pydamage",
                "VN": __version__,
                "CL": " ".join(sys.argv),
            }
        )
        refs = al.references
        with pysam.AlignmentFile(outname, "wb", threads=threads, header=hd) as out:
            for ref in tqdm(refs, desc="Rescaling quality scores"):
                if grouped:
                    pydam_ref = "reference"
                else:
                    pydam_ref = ref
                dmg = damage_model()
                if pydam_ref in read_dict:
                    if (
                        threshold
                        and 
                        threshold <= damage_dict["predicted_accuracy"][pydam_ref]
                    ) and (alpha and alpha >= damage_dict["qvalue"][pydam_ref]):
                        dmg_pmf = dmg.fit(
                            x=np.arange(wlen),
                            p=damage_dict["damage_model_p"][pydam_ref],
                            pmin=damage_dict["damage_model_pmin"][pydam_ref],
                            pmax=damage_dict["damage_model_pmax"][pydam_ref],
                        )
                        for read in al.fetch(ref):
                            if read.query_name in read_dict[pydam_ref]:
                                qual = read.query_qualities
                                read.query_sequence = read.query_sequence
                                read.query_qualities = array(
                                    "B",
                                    rescale_qual(
                                        qual,
                                        dmg_pmf,
                                        read_dict[pydam_ref][read.query_name],
                                        reverse=read.is_reverse,
                                    ),
                                )
                            out.write(read)
                    else:
                        for read in al.fetch(ref):
                            out.write(read)
                else:
                    for read in al.fetch(ref):
                        out.write(read)
