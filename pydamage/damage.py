#!/usr/bin/env python3

import pysam
from pydamage.parse_damage import damage_al
from pydamage.model_fit import fit_models
from pydamage import models
import numpy as np
from numba import jit
from tqdm import tqdm

def sort_count_array_dict(int_array):
    """Sorts and counts unique values in an array

    Args:
        int_array (np.array): Array of integers
    Returns:
        np.array: sorted unique values
        np.array: counts of unique values
    """
    pos, counts = np.unique(int_array, return_counts=True)
    return pos, counts

class al_to_damage:
    def __init__(self, reference, al_handle, wlen, g2a):
        """Constructor of the class

        Args:
            reference (string): a reference from the indexed bam file
            al_handle(pysam.AlignmentFile)
            wlen (int): window length
            g2a (bool): use GtoA transitions
        """
        self.alignments = al_handle.fetch(reference)
        self.reference = reference
        self.wlen = wlen
        self.g2a = g2a
        # self.alignments = al_file

    # def __repr__(self):
    #     return(f"Reference {reference}")

    def get_damage(self, show_al):
        """Compute CtoT substitutions

        Args:
            show_al(bool): print alignments representations
        Returns:
            list: C - positions from 5' of C in reads
            list: CT- positions from 5' of CtoT transitions in reads
            list: G - positions from 3' of G in reverse reads
            list: GA- positions from 3' of GtoA transitions in reverse reads
            list: damage_bases - positions of C2T and G2A bases (if reversed)
            list: C_G_bases - positions of C2C and G2G bases (if reversed)
        """
        self.C = []
        self.CT = []
        self.damage_bases = []
        self.C_G_bases = []
        self.G = []
        self.GA = []
        self.no_mut = []
        self.read_dict = {self.reference: dict()}
        if self.reference is None:
            iterator = tqdm(self.alignments)
        else:
            iterator = self.alignments
        for al in iterator:
            if al.is_unmapped is False:
                all_damage = damage_al(
                    reference=al.get_reference_sequence(),
                    read_name=al.query_name,
                    is_reverse=al.is_reverse,
                    ref_name=self.reference,
                    query=al.query_sequence,
                    cigartuple=al.cigartuples,
                    wlen=self.wlen,
                    g2a=self.g2a,
                    show_al=show_al,
                )
                self.C += all_damage["C"]
                self.CT += all_damage["CT"]
                self.G += all_damage["G"]
                self.GA += all_damage["GA"]
                self.C_G_bases += all_damage["C"]
                CT_GA = all_damage["CT"]
                if self.g2a:
                    CT_GA += all_damage["GA"]
                    self.C_G_bases += all_damage["G"]
                self.damage_bases += CT_GA
                self.no_mut += all_damage["no_mut"]
                if len(CT_GA) > 0 and (
                    (al.is_reverse and self.g2a) or (not al.is_reverse)
                ):
                    self.read_dict[self.reference][al.query_name] = np.array(CT_GA)

    def compute_damage(self):
        """Computes the amount of damage for statistical modelling"""

        # All C in reference
        C_pos, C_counts = sort_count_array_dict(np.array(self.C, dtype="uint32"))

        # All G in reference
        G_pos, G_counts = sort_count_array_dict(np.array(self.G, dtype="uint32"))

        # CtoT transitions
        CT_pos, CT_counts = sort_count_array_dict(np.array(self.CT, dtype="uint32"))

        # GtoA transitions
        GA_pos, GA_counts = sort_count_array_dict(np.array(self.GA, dtype="uint32"))

        # All transitions
        damage_bases_pos, damage_bases_counts = sort_count_array_dict(np.array(self.damage_bases, dtype="uint32"))

        # All C and G in reference
        C_G_bases_pos, C_G_bases_counts = sort_count_array_dict(np.array(self.C_G_bases, dtype="uint32"))

        # All conserved C and G
        no_mut_pos, no_mut_counts = sort_count_array_dict(np.array(self.no_mut, dtype="uint32"))

        CT_damage_amount = np.zeros(self.wlen)
        CT_damage_amount[CT_pos] = CT_counts / C_counts[CT_pos]

        GA_damage_amount = np.zeros(self.wlen)
        GA_damage_amount[GA_pos] = GA_counts / G_counts[GA_pos]

        damage_amount = np.zeros(self.wlen)
        damage_amount[damage_bases_pos] = damage_bases_counts / C_G_bases_counts[damage_bases_pos]

        _ = np.zeros(self.wlen, dtype="uint32")
        _[damage_bases_pos] = damage_bases_counts
        damage_bases_counts = _

        _ = np.zeros(self.wlen, dtype="uint32")
        _[no_mut_pos] = no_mut_counts
        no_mut_counts = _

        return (
            damage_bases_counts,  # Number of CtoT and GtoA per position
            no_mut_counts,  # Number of conserved C and G per position
            CT_damage_amount,
            GA_damage_amount,
            damage_amount,
        )


@jit
def avg_coverage_contig(pysam_cov):
    """Computes average coverage of a reference

    Args:
        pysam_cov (np.array): Four dimensional array of coverage for each base
    Returns:
        float: mean coverage of reference
    """
    return (
        np.mean(
            np.sum(pysam_cov, axis=0)
        )
    )


def check_model_fit(model_dict, wlen, verbose):
    """Check if model fitting makes sense

    Args:
        model_dict (dict): Dictionary containing LR test results
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

    return model_dict


def test_damage(ref, bam, mode, wlen, g2a, show_al, process, verbose):
    """Prepare data and run LRtest to test for damage

    Args:
        ref (str): name of reference in alignment file
        bam (str): bam file
        mode (str): opening mode of alignment file
        wlen (int): window length
        g2a (bool): Use GtoA transitions
        show_al (bool): Show alignment representations
        process (int): Number of process for parallelization
        verbose (bool): Run in verbose mode
    Returns:
        dict: Dictionary containing LR test results
    """
    al_handle = pysam.AlignmentFile(bam, mode=mode, threads=process)
    try:
        if ref is None:
            all_references = al_handle.references
            cov = avg_coverage_contig(
                np.concatenate(
                    [np.array(al_handle.count_coverage(contig=ref), dtype="uint16") for ref in all_references], 
                    axis=1
                )
            )
            nb_reads_aligned = np.sum(
                [al_handle.count(contig=ref) for ref in all_references]
            )
            reflen = np.sum(
                [al_handle.get_reference_length(ref) for ref in all_references]
            )
            refname = "reference"
        else:
            cov = avg_coverage_contig(np.array(al_handle.count_coverage(contig=ref), dtype="uint16"))
            nb_reads_aligned = al_handle.count(contig=ref)
            reflen = al_handle.get_reference_length(ref)
            refname = ref

        al = al_to_damage(reference=ref, al_handle=al_handle, wlen=wlen, g2a=g2a)
        al.get_damage(show_al=show_al)
        read_dict = al.read_dict
        if len(read_dict.keys()) == 1 and list(read_dict.keys())[0] is None:
            read_dict[refname] = read_dict.pop(None)
        (
            mut_count,
            conserved_count,
            CT_damage,
            GA_damage,
            all_damage,
        ) = al.compute_damage()

        model_A = models.damage_model()
        model_B = models.null_model()
        test_res = fit_models(
            ref=ref,
            model_A=model_A,
            model_B=model_B,
            damage=all_damage,
            mut_count=mut_count,
            conserved_count=conserved_count,
            verbose=verbose,
        )
        test_res["reference"] = refname
        test_res["nb_reads_aligned"] = nb_reads_aligned
        test_res["coverage"] = cov
        test_res["reflen"] = reflen

        CT_log = {}
        GA_log = {}

        for i in range(wlen):
            CT_log[f"CtoT-{i}"] = CT_damage[i]
            GA_log[f"GtoA-{i}"] = GA_damage[i]
        test_res.update(CT_log)
        test_res.update(GA_log)

        return (check_model_fit(test_res, wlen, verbose), read_dict)

    except (ValueError, RuntimeError) as e:
        if verbose:
            print(f"Model fitting for {ref} failed")
            print(f"Model fitting error: {e}")
            print(
                f"nb_reads_aligned: {nb_reads_aligned} - coverage: {cov}"
                f" - reflen: {reflen}\n"
            )
        return False


# def get_damage_group(get_damage_list):
#     """Prepare data and run LR test to test for damage on reference grouped
#     as one.

#     Args:
#         get_damage_list (list): results of get_damage as list of dict
#     Returns:
#         list:[ct_data, ga_data, cc_data, c_data, g_data, all_bases]
#     """

#     al = al_to_damage(reference=ref, al_handle=al_handle)
#     (
#         CT,
#         GA,
#         damage_bases,
#         forward,
#         reverse,
#     ) = al.get_damage(wlen=wlen, show_al=show_al)

#     return [
#         CT,
#         GA,
#         damage_bases,
#         forward,
#         reverse,
#         cov,
#         nb_reads_aligned,
#         reflen,
#     ]


# def test_damage_group(
#     ct_data,
#     ga_data,
#     damage_bases,
#     forward,
#     reverse,
#     nb_reads_aligned,
#     cov,
#     reflen,
#     wlen,
#     verbose,
# ):
#     """Performs damage test

#     Args:
#         ct_data (list of int): List of positions from 5' with CtoT transitions
#         ga_data (list of int): List of positions from 3' with GtoA transitions
#         damage_bases (list of int): List of positions with damage
#         forward (list of int): List of positions from 5' where a base is aligned
#         reverse(list of int): List of positions from 3' where a base is aligned
#         nb_reads_aligned(int): number of reads aligned
#         cov(float): average coverage across all references
#         reflen(int): length of all references
#         wlen (int): window length
#         verbose(bool): Verbose
#     Returns:
#         dict: Dictionary containing LR test results
#     """
#     ref = "reference"
#     try:
#         if ct_data:
#             model_A = models.damage_model()
#             model_B = models.null_model()
#             test_res = fit_models(
#                 ref="reference",
#                 model_A=model_A,
#                 model_B=model_B,
#                 ct_data=ct_data,
#                 cc_data=cc_data,
#                 ga_data=ga_data,
#                 all_bases=all_bases,
#                 wlen=wlen,
#                 verbose=verbose,
#             )
#             test_res["reference"] = ref
#             test_res["nb_reads_aligned"] = nb_reads_aligned
#             test_res["coverage"] = cov
#             test_res["reflen"] = reflen

#             return check_model_fit(test_res, wlen, verbose)

#     except (ValueError, RuntimeError) as e:
#         if verbose:
#             print(f"Model fitting for {ref} failed")
#             print(f"Model fitting error: {e}")
#             print(
#                 f"nb_reads_aligned: {nb_reads_aligned} - coverage: {cov} "
#                 "- reflen: {reflen}\n"
#             )
#         return False
