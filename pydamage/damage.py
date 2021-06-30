#!/usr/bin/env python3

import pysam
from pydamage.parse_damage import damage_al
from pydamage.model_fit import fit_models
from pydamage import models
import numpy as np


class al_to_damage:
    def __init__(self, reference, al_handle, wlen):
        """Constructor of the class

        Args:
            reference (string): a reference from the indexed bam file
            al_handle(pysam.AlignmentFile)
            wlen (int): window length

        """
        self.alignments = al_handle.fetch(reference)
        self.reference = reference
        self.wlen = wlen
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
        self.no_mut = []
        for al in self.alignments:
            if al.is_unmapped is False:
                all_damage = damage_al(
                    reference=al.get_reference_sequence(),
                    read_name=al.query_name,
                    is_reverse=al.is_reverse,
                    ref_name=self.reference,
                    query=al.query_sequence,
                    cigartuple=al.cigartuples,
                    wlen=self.wlen,
                    show_al=show_al,
                )
                self.C += all_damage["C"]
                self.CT += all_damage["CT"]
                self.damage_bases += all_damage["CT"]
                self.C_G_bases += all_damage["C"]
                self.no_mut += all_damage["no_mut"]

    def compute_damage(self):
        """Computes the amount of damage for statistical modelling"""

        # All C in reference
        C_pos, C_counts = np.unique(
            np.sort(
                self.C,
            ),
            return_counts=True,
        )
        C_dict = dict(zip(C_pos, C_counts))

        # CtoT transitions
        CT_pos, CT_counts = np.unique(np.sort(self.CT), return_counts=True)
        CT_dict = dict(zip(CT_pos, CT_counts))

        # All transitions
        damage_bases_pos, damage_bases_counts = np.unique(
            np.sort(self.damage_bases), return_counts=True
        )
        damage_bases_dict = dict(zip(damage_bases_pos, damage_bases_counts))

        # All C and G in reference
        C_G_bases_pos, C_G_bases_counts = np.unique(
            np.sort(self.C_G_bases), return_counts=True
        )
        C_G_bases_dict = dict(zip(C_G_bases_pos, C_G_bases_counts))

        # All conserved C and G
        no_mut_pos, no_mut_counts = np.unique(np.sort(self.no_mut), return_counts=True)
        no_mut_dict = dict(zip(no_mut_pos, no_mut_counts))

        CT_damage_amount = []
        GA_damage_amount = []
        damage_amount = []

        # Checking for non covered positions
        for i in range(self.wlen):
            if i not in CT_dict:
                CT_dict[i] = 0
            if i not in damage_bases_dict:
                damage_bases_dict[i] = 0
            if i not in no_mut_dict:
                no_mut_dict[i] = 0
            if i not in C_dict:
                CT_damage_amount.append(0)
            else:
                CT_damage_amount.append(CT_dict[i] / C_dict[i])
            if i not in C_G_bases_dict:
                damage_amount.append(0)
            else:
                damage_amount.append(damage_bases_dict[i] / C_G_bases_dict[i])

        return (
            np.array(
                list(damage_bases_dict.values())
            ),  # Number of CtoT and GtoA per position
            np.array(
                list(no_mut_dict.values())
            ),  # Number of conserved C and G per position
            np.array(CT_damage_amount),
            np.array(GA_damage_amount),
            np.array(damage_amount),
        )


def avg_coverage_contig(pysam_cov):
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

    # Check that all the first wlen bases are covered
    # if len(model_dict["base_cov"]) < wlen:
    #     return False
    # if np.any(np.array(model_dict["base_cov"][:wlen]) == 0):
    #     if verbose:
    #         print(
    #             f"Could not reliably fit a model to {model_dict['reference']}"
    #             "because of too few reads aligned"
    #         )
    #     return False
    return model_dict


def test_damage(ref, bam, mode, wlen, show_al, process, verbose):
    """Prepare data and run LRtest to test for damage

    Args:
        ref (str): name of referene in alignment file
        bam (str): bam file
        mode (str): opening mode of alignment file
        wlen (int): window length
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
            cov = np.mean(
                [
                    avg_coverage_contig(al_handle.count_coverage(contig=ref))
                    for ref in all_references
                ]
            )
            nb_reads_aligned = np.sum(
                [al_handle.count(contig=ref) for ref in all_references]
            )
            reflen = np.sum(
                [al_handle.get_reference_length(ref) for ref in all_references]
            )
            refname = "reference"
        else:
            cov = avg_coverage_contig(al_handle.count_coverage(contig=ref))
            nb_reads_aligned = al_handle.count(contig=ref)
            reflen = al_handle.get_reference_length(ref)
            refname = ref

        al = al_to_damage(reference=ref, al_handle=al_handle, wlen=wlen)
        al.get_damage(show_al=show_al)
        (
            mut_count,
            conserved_count,
            CT_damage,
            GA_damage,
            all_damage,
        ) = al.compute_damage()
        # if all_damage:
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
        test_res.update(CT_log)
        test_res.update(GA_log)

        # for i in range(qlen):
        #     if i not in ydata_counts:
        #         ydata_counts[i] = np.nan
        #     if f"CtoT-{i}" not in ctot_out:
        #         ctot_out[f"CtoT-{i}"] = np.nan
        #     if f"GtoA-{i}" not in gtoa_out:
        #         gtoa_out[f"GtoA-{i}"] = np.nan

        # print(test_res)

        return check_model_fit(test_res, wlen, verbose)

    except (ValueError, RuntimeError) as e:
        if verbose:
            print(f"Model fitting for {ref} failed")
            print(f"Model fitting error: {e}")
            print(
                f"nb_reads_aligned: {nb_reads_aligned} - coverage: {cov}"
                " - reflen: {reflen}\n"
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
