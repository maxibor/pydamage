#!/usr/bin/env python3

from array import array
import sys
import numpy as np
import pysam
from tqdm import tqdm
from numba import jit

from pydamage.models import damage_model
from pydamage import __version__


@jit(cache=True, fastmath=True, nopython=True)
def phred_to_prob(qual: np.ndarray) -> np.ndarray:
    """Convert Phred quality score to probability."""
    return 10.0 ** (-qual / 10.0)


@jit(cache=True, fastmath=True, nopython=True)
def compute_new_prob(e: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Compute new probability of base calling error accounting for ancient damage."""
    return np.round(-10.0 * np.log10(1.0 - (1.0 - e) * (1.0 - d))).astype(np.int64)


def rescale_qual(read_qual: list, dmg_pmf: np.ndarray, damage_bases: list, reverse: bool) -> np.ndarray:
    """Rescale quality scores using damage model PMF."""
    e = phred_to_prob(np.asarray(read_qual, dtype=np.float64))
    if reverse:
        e = e[::-1]
    d = np.zeros(len(read_qual), dtype=np.float64)
    dmg_arr = np.asarray(damage_bases, dtype=np.int64)
    if len(dmg_arr) > 0:
        valid_dmg = dmg_arr[dmg_arr < len(dmg_pmf)]
        valid_dmg = valid_dmg[valid_dmg < len(d)]
        d[valid_dmg] = dmg_pmf[valid_dmg]
    r = compute_new_prob(e, d)
    if reverse:
        r = r[::-1]
    return r


def rescale_bam(
    bam: str,
    threshold: float,
    alpha: float,
    wlen: int,
    damage_dict: dict,
    read_dict: dict,
    grouped: bool,
    outname: str,
    threads: int = 1,
):
    """Rescale quality scores in BAM file using damage model parameters."""
    with pysam.AlignmentFile(bam, "rb", threads=threads) as al:
        hd = al.header.to_dict()
        if "PG" not in hd:
            hd["PG"] = []
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
                pydam_ref = "reference" if grouped else ref
                dmg = damage_model()
                if pydam_ref in read_dict:
                    pred_acc = damage_dict.get("predicted_accuracy", {}).get(pydam_ref, 0.0)
                    qval = damage_dict.get("qvalue", {}).get(pydam_ref, 1.0)
                    if (threshold is None or threshold <= pred_acc) and (alpha is None or alpha >= qval):
                        dmg_pmf = dmg.fit(
                            x=np.arange(wlen),
                            p=damage_dict["damage_model_p"][pydam_ref],
                            pmin=damage_dict["damage_model_pmin"][pydam_ref],
                            pmax=damage_dict["damage_model_pmax"][pydam_ref],
                            wlen=wlen,
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
