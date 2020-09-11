#!/usr/bin/env python3

import pysam
from pydamage.parse_damage import damage_al
from pydamage.model_fit import fit_models
from pydamage import models
from pydamage.gc_content import fastaFile
import numpy as np


class al_to_damage:
    def __init__(self, reference, al_handle):
        """Constructor of the class

        Args:
            reference (string): a reference from the indexed bam file
            al_handle(pysam.AlignmentFile)
        """
        self.alignments = al_handle.fetch(reference)

        # self.alignments = al_file

    # def __repr__(self):
    #     return(f"Reference {reference}")

    def get_damage(self, wlen, show_al):
        """Compute CtoT substitutions

        Args:
            wlen (int): window length
            show_al(bool): print alignments representations
        Returns:
            list: all_ct - positions of CtoT transitions in reads
            list: all_ga - positions of GtoA transitions in reads
            list: all_cc - positions of C in ref and reads, in reads
            list: all_c - positions of C in reads
            list: all_g - positions of G in reads
            list: all_bases - position of all bases
        """

        all_ct = []
        all_ga = []
        all_cc = []
        all_c = []
        all_g = []
        all_bases = []
        for al in self.alignments:
            if al.is_reverse == False and al.is_unmapped == False:
                cigar = al.cigartuples
                ref = al.get_reference_sequence()
                quer = al.query_sequence

                all_damage = damage_al(
                    reference=ref,
                    query=quer,
                    cigartuple=cigar,
                    wlen=wlen,
                    show_al=show_al,
                )
                all_ct += all_damage["CT"]
                all_ga += all_damage["GA"]
                all_cc += all_damage["CC"]
                all_c += all_damage["C"]
                all_g += all_damage["G"]
                all_bases += all_damage["all"]

        return (all_ct, all_ga, all_cc, all_c, all_g, all_bases)


def avg_coverage(pysam_cov):
    """Computes average coverage of a reference

    Args:
        pysam_cov (np.array): Four dimensional array of coverage for each base
    Returns:
        float: mean coverage of reference
    """
    A = np.array(pysam_cov[0], dtype=int)
    C = np.array(pysam_cov[1], dtype=int)
    G = np.array(pysam_cov[2], dtype=int)
    T = np.array(pysam_cov[3], dtype=int)
    cov_all_bases = A + C + G + T
    cov = np.mean(cov_all_bases)
    return cov


def check_model_fit(model_dict, wlen, verbose):
    """Check if model fitting makes sense

    Args:
        model_dict (dict): Dictionary containing Vuong test results
        wlen (int): window length
        verbose (bool): Run in verbose mode
    Returns:
        bool or model_dict(dict): False if test fails, model_dict otherwise
    """
    # Check that no fitted parameters or stdev are infinite
    if np.isinf(np.array(model_dict["model_params"])).any():
        if verbose:
            print(f"Unreliable model fit for {model_dict['reference']}")
        return False

    # Check that all the first wlen bases are covered
    if np.any(np.array(model_dict["base_cov"][:wlen]) == 0):
        if verbose:
            print(
                f"Could not reliably fit a model to {model_dict['reference']} because of too few reads aligned"
            )
        return False
    return model_dict


def test_damage(
    ref, bam, fasta, mode, wlen, show_al, min_al, min_cov, process, verbose
):
    """Prepare data and run Vuong's test to test for damage

    Args:
        ref (str): name of referene in alignment file
        bam (str): bam file
        fasta(str): Path to fasta file containing reference sequences
        mode (str): opening mode of alignment file
        wlen (int): window length
        show_al (bool): Show alignment representations
        min_al (int): Minimum  number of aligned reads
        min_cov (float): Minimum coverage
        process (int): Number of process for parallelization
        verbose (bool): Run in verbose mode
    Returns:
        dict: Dictionary containing Vuong test results
    """
    al_handle = pysam.AlignmentFile(bam, mode=mode, threads=process)
    fasta_handle = fastaFile(fasta)
    fasta_gc = fasta_handle.compute_gc()
    try:
        cov = avg_coverage(al_handle.count_coverage(contig=ref))
        nb_reads_aligned = al_handle.count(contig=ref)
        reflen = al_handle.get_reference_length(ref)

        if nb_reads_aligned >= min_al or cov >= min_cov:
            al = al_to_damage(reference=ref, al_handle=al_handle)
            ct_data, ga_data, cc_data, c_data, g_data, all_bases = al.get_damage(
                wlen=wlen, show_al=show_al
            )
            if ct_data:
                model_A = models.damage_model()
                model_B = models.null_model()
                test_res = fit_models(
                    ref=ref,
                    model_A=model_A,
                    model_B=model_B,
                    ct_data=ct_data,
                    cc_data=cc_data,
                    ga_data=ga_data,
                    all_bases=all_bases,
                    wlen=wlen,
                    verbose=verbose,
                )
                test_res["reference"] = ref
                test_res["nb_reads_aligned"] = nb_reads_aligned
                test_res["coverage"] = cov
                test_res["reflen"] = reflen
                test_res["gc_content"] = fasta_gc[ref]

                return check_model_fit(test_res, wlen, verbose)
        else:
            if verbose:
                print(
                    f"Did not attempt to fit a model to {ref} because of too few reads aligned"
                )
                print(
                    f"nb_reads_aligned: {nb_reads_aligned} - coverage: {cov} - reflen: {reflen}\n"
                )
            pass
    except (ValueError, RuntimeError) as e:
        if verbose:
            print(f"Model fitting for {ref} failed")
            print(f"Model fitting error: {e}")
            print(
                f"nb_reads_aligned: {nb_reads_aligned} - coverage: {cov} - reflen: {reflen}\n"
            )
        return False
