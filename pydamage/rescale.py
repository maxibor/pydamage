import pysam
import numpy as np
from array import array
from pydamage.models import damage_model
from tqdm import tqdm
from numba import jit


@jit
def phred_to_prob(qual):
    """Convert Phred quality score to probability

    Args:
        qual (array): Array of unsigned integer Phred quality scores

    Returns:
        np.array(int): Array of read error probabilities
    """
    return 10 ** (-qual / 10)


@jit
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


def rescale_bam(bam, threshold, alpha, damage_dict, read_dict, outname):
    """Rescale quality scores in BAM file using damage model

    Args:
        bam (str): Path to BAM file
        threshold (float): Predicted accuracy threshold
        alpha (float): Q-value threshold
        damage_dict (dict): Damage model parameters
        read_dict (dict): Dictionary of read names
        outname (str): Path to output BAM file
    """
    with pysam.AlignmentFile(bam, "rb") as al:
        refs = al.references
        with pysam.AlignmentFile(outname, "wb", template=al) as out:
            for ref in tqdm(refs, desc="Rescaling quality scores"):
                dmg = damage_model()
                if ref in read_dict:
                    pass_filter = False
                    if (
                        threshold
                        and threshold <= damage_dict["predicted_accuracy"][ref]
                    ) and (alpha and alpha >= damage_dict["qvalue"][ref]):
                        pass_filter = True
                    if pass_filter:
                        dmg_pmf = dmg.fit(
                            x=np.arange(400),
                            p=damage_dict["damage_model_p"][ref],
                            pmin=damage_dict["damage_model_pmin"][ref],
                            pmax=damage_dict["damage_model_pmax"][ref],
                        )
                        for read in al.fetch(ref):
                            if read.query_name in read_dict[ref]:
                                qual = read.query_qualities
                                read.query_qualities = array(
                                    "B",
                                    rescale_qual(
                                        qual,
                                        dmg_pmf,
                                        read_dict[ref][read.query_name],
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
