import pysam
import numpy as np
from array import array
from pydamage.models import damage_model


def phred_to_prob(qual):
    """Convert Phred quality score to probability

    Args:
        qual (array): Array of unsigned integer Phred quality scores

    Returns:
        np.array(int): Array of read error probabilities
    """
    return 10 ** (-np.array(qual) / 10)


def rescale_qual(read_qual, dmg_pmf, damage_bases):
    """Rescale quality scores using damage model
    Args:
        read_qual (array): Array of Phred quality scores
        dmg_pmf (array): Array of damage model probabilities
        damage_bases (array): Array of positions with damage
    Returns:
        np.array(int): Array of rescaled Phred quality scores
    """
    e = phred_to_prob(read_qual)
    d = np.zeros(len(read_qual))
    d[damage_bases] = dmg_pmf[damage_bases]
    return array(
        "B", np.round(-10 * np.log10(1 - np.multiply(1 - e, 1 - d)), 0).astype(int)
    )


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
    damage_dict = {v["reference"]: v for v in damage_dict}
    print(damage_dict)
    with pysam.AlignmentFile(bam, "rb") as al:
        refs = al.references
        with pysam.AlignmentFile(outname, "wb", template=al) as out:
            for ref in refs:
                dmg = damage_model()
                if ref in read_dict:
                    pass_filter = False
                    if (
                        threshold
                        and threshold <= damage_dict[ref]["predicted_accuracy"]
                    ) and (alpha and alpha >= damage_dict[ref]["qvalue"]):
                        pass_filter = True
                    if pass_filter:
                        dmg_pmf = dmg.fit(x=np.arange(400), **damage_dict[ref])
                        for read in al.fetch(ref):
                            if read.query_name in read_dict[ref]:
                                qual = read.query_qualities
                                read.query_qualities = rescale_qual(
                                    qual, dmg_pmf, read_dict[ref][read.query_name]
                                )
                                print(
                                    np.mean(np.array(qual)),
                                    np.mean(
                                        np.array(
                                            rescale_qual(
                                                qual,
                                                dmg_pmf,
                                                read_dict[ref][read.query_name],
                                            )
                                        )
                                    ),
                                )

                            out.write(read)
                    else:
                        for read in al.fetch(ref):
                            out.write(read)
                else:
                    for read in al.fetch(ref):
                        out.write(read)
