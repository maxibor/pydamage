#!/usr/bin/env python3

import logging
import numpy as np
import pysam
from tqdm import tqdm

from pydamage.parse_damage import damage_al
from pydamage.model_fit import fit_models
from pydamage.contig_stats import get_contig_stats, get_ref_stats
from pydamage import models

# Thread-local/process-local cached BAM handle
_process_local_al_handle = None


def init_worker_handle(bam: str, mode: str):
    """Initializer for multiprocessing worker processes to keep an open BAM handle."""
    global _process_local_al_handle
    _process_local_al_handle = pysam.AlignmentFile(bam, mode=mode, threads=1)


def close_worker_handle():
    """Closes worker persistent BAM handle."""
    global _process_local_al_handle
    if _process_local_al_handle is not None:
        try:
            _process_local_al_handle.close()
        except Exception:
            pass
        _process_local_al_handle = None


def sort_count_array_dict(int_array: np.ndarray) -> tuple:
    """Sorts and counts unique values in an array.

    Args:
        int_array (np.ndarray): Array of integers
    Returns:
        np.ndarray: sorted unique values
        np.ndarray: counts of unique values
    """
    pos, counts = np.unique(int_array, return_counts=True)
    return pos, counts


class al_to_damage:
    """Extracts and computes damage substitution patterns and stats from alignments in a single pass."""

    def __init__(self, reference, al_handle: pysam.AlignmentFile, wlen: int, g2a: bool, subsample):
        self.alignments = al_handle.fetch(reference)
        self.reference = reference
        self.al_handle = al_handle
        self.wlen = wlen
        self.g2a = g2a
        self.subsample = subsample

        self.C = []
        self.CT = []
        self.damage_bases = []
        self.C_G_bases = []
        self.G = []
        self.GA = []
        self.no_mut = []
        self.read_dict = {self.reference: {}}

        self.nb_reads_aligned = 0
        self.cov_sum = 0
        self.reflen = 0

    def __repr__(self):
        return f"Reference {self.reference}"

    def get_damage(self, show_al: bool = False, need_read_dict: bool = True):
        """Extract substitution positions and calculate alignment coverage in a single pass."""
        self.C = []
        self.CT = []
        self.damage_bases = []
        self.C_G_bases = []
        self.G = []
        self.GA = []
        self.no_mut = []
        self.read_dict = {self.reference: {}} if need_read_dict else {}
        self.nb_reads_aligned = 0
        self.cov_sum = 0

        if self.reference is not None and hasattr(self.al_handle, "get_reference_length"):
            try:
                self.reflen = self.al_handle.get_reference_length(self.reference)
            except Exception:
                self.reflen = 0

        iterator = (
            tqdm(self.alignments, desc="Compute damage for entire reference")
            if self.reference is None
            else self.alignments
        )

        for al in iterator:
            if al.is_unmapped:
                continue

            self.nb_reads_aligned += 1
            cigartuples = al.cigartuples
            quals = al.query_qualities

            # Compute quality-filtered coverage sum in the same pass
            if cigartuples and quals is not None:
                qpos = 0
                for op, length in cigartuples:
                    if op in (0, 7, 8):
                        # bases with quality >= 15 match pysam count_coverage default
                        self.cov_sum += sum(1 for q in quals[qpos : qpos + length] if q >= 15)
                        qpos += length
                    elif op in (1, 4):
                        qpos += length

            if not self.subsample or np.random.rand() < self.subsample:
                all_damage = damage_al(
                    reference=al.get_reference_sequence(),
                    read_name=al.query_name,
                    is_reverse=al.is_reverse,
                    ref_name=self.reference,
                    query=al.query_sequence,
                    cigartuple=cigartuples if cigartuples else [],
                    wlen=self.wlen,
                    g2a=self.g2a,
                    show_al=show_al,
                )

                c_list = all_damage["C"]
                ct_list = all_damage["CT"]
                g_list = all_damage["G"]
                ga_list = all_damage["GA"]
                no_mut_list = all_damage["no_mut"]

                self.C.extend(c_list)
                self.CT.extend(ct_list)
                self.G.extend(g_list)
                self.GA.extend(ga_list)
                self.C_G_bases.extend(c_list)

                ct_ga = list(ct_list)
                if self.g2a:
                    ct_ga.extend(ga_list)
                    self.C_G_bases.extend(g_list)

                self.damage_bases.extend(ct_ga)
                self.no_mut.extend(no_mut_list)

                if need_read_dict and ct_ga and ((al.is_reverse and self.g2a) or not al.is_reverse):
                    self.read_dict[self.reference][al.query_name] = np.array(ct_ga)

    def compute_damage(self) -> tuple:
        """Computes damage frequencies and counts per position across the window."""
        damage_bases_counts = np.zeros(self.wlen, dtype="uint32")
        no_mut_counts = np.zeros(self.wlen, dtype="uint32")
        CT_damage_amount = np.zeros(self.wlen, dtype=np.float64)
        GA_damage_amount = np.zeros(self.wlen, dtype=np.float64)
        damage_amount = np.zeros(self.wlen, dtype=np.float64)

        if self.C:
            C_pos, C_counts = sort_count_array_dict(np.array(self.C, dtype="uint32"))
            if self.CT:
                CT_pos, CT_counts = sort_count_array_dict(np.array(self.CT, dtype="uint32"))
                C_pos_dict = {pos: i for i, pos in enumerate(C_pos)}
                CT_pos_indices = [C_pos_dict[pos] for pos in CT_pos if pos in C_pos_dict]
                valid_ct_pos = [pos for pos in CT_pos if pos in C_pos_dict]
                CT_damage_amount[valid_ct_pos] = CT_counts[:len(valid_ct_pos)] / C_counts[CT_pos_indices]

        if self.G:
            G_pos, G_counts = sort_count_array_dict(np.array(self.G, dtype="uint32"))
            if self.GA:
                GA_pos, GA_counts = sort_count_array_dict(np.array(self.GA, dtype="uint32"))
                G_pos_dict = {pos: i for i, pos in enumerate(G_pos)}
                GA_pos_indices = [G_pos_dict[pos] for pos in GA_pos if pos in G_pos_dict]
                valid_ga_pos = [pos for pos in GA_pos if pos in G_pos_dict]
                GA_damage_amount[valid_ga_pos] = GA_counts[:len(valid_ga_pos)] / G_counts[GA_pos_indices]

        if self.C_G_bases and self.damage_bases:
            C_G_bases_pos, C_G_bases_counts = sort_count_array_dict(np.array(self.C_G_bases, dtype="uint32"))
            damage_bases_pos, damage_bases_counts_raw = sort_count_array_dict(np.array(self.damage_bases, dtype="uint32"))
            CG_pos_dict = {pos: i for i, pos in enumerate(C_G_bases_pos)}
            CG_pos_indices = [CG_pos_dict[pos] for pos in damage_bases_pos if pos in CG_pos_dict]
            valid_dmg_pos = [pos for pos in damage_bases_pos if pos in CG_pos_dict]
            damage_amount[valid_dmg_pos] = damage_bases_counts_raw[:len(valid_dmg_pos)] / C_G_bases_counts[CG_pos_indices]

            for pos, count in zip(damage_bases_pos, damage_bases_counts_raw):
                if pos < self.wlen:
                    damage_bases_counts[pos] = count

        if self.no_mut:
            no_mut_pos, no_mut_counts_raw = sort_count_array_dict(np.array(self.no_mut, dtype="uint32"))
            for pos, count in zip(no_mut_pos, no_mut_counts_raw):
                if pos < self.wlen:
                    no_mut_counts[pos] = count

        return (
            damage_bases_counts,
            no_mut_counts,
            CT_damage_amount,
            GA_damage_amount,
            damage_amount,
        )


def check_model_fit(model_dict: dict, wlen: int, verbose: bool = False):
    """Check if model fitting parameters are valid and finite."""
    if not isinstance(model_dict, dict):
        return False
    if np.isinf(np.array(model_dict["model_params"])).any():
        if verbose:
            logging.warning(f"Unreliable model fit for {model_dict.get('reference')}")
        return False
    return model_dict


def test_damage_with_handle(
    ref: str,
    al_handle: pysam.AlignmentFile,
    wlen: int,
    g2a: bool,
    subsample: float,
    show_al: bool,
    verbose: bool,
    need_read_dict: bool = True,
) -> tuple:
    """Core damage testing function reusing an open AlignmentFile handle."""
    reflen, cov, nb_reads_aligned = 0, 0.0, 0
    read_dict = {}

    try:
        if ref is None:
            if verbose:
                logging.info("Computing alignment stats for entire reference")
            reflen, cov, nb_reads_aligned = get_ref_stats(al_handle.filename.decode() if isinstance(al_handle.filename, bytes) else str(al_handle.filename))
            refname = "reference"
        else:
            refname = ref

        al = al_to_damage(
            reference=ref, al_handle=al_handle, wlen=wlen, g2a=g2a, subsample=subsample
        )
        al.get_damage(show_al=show_al, need_read_dict=need_read_dict)
        read_dict = al.read_dict

        if ref is not None:
            reflen = al.reflen
            nb_reads_aligned = al.nb_reads_aligned
            cov = al.cov_sum / reflen if reflen > 0 else 0.0

        if subsample:
            cov = cov * subsample
            nb_reads_aligned = int(np.round(nb_reads_aligned * subsample, 0))

        if len(read_dict) == 1 and None in read_dict:
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

        for i in range(wlen):
            test_res[f"CtoT-{i}"] = CT_damage[i]
            test_res[f"GtoA-{i}"] = GA_damage[i]

        return check_model_fit(test_res, wlen, verbose), read_dict

    except (ValueError, RuntimeError) as e:
        if verbose:
            logging.warning(f"Model fitting for {ref} failed: {e}")
            logging.warning(
                f"nb_reads_aligned: {nb_reads_aligned} - coverage: {cov} - reflen: {reflen}\n"
            )
        return False, read_dict


def test_damage_worker_task(
    ref: str,
    wlen: int,
    g2a: bool,
    subsample: float,
    show_al: bool,
    verbose: bool,
    need_read_dict: bool,
) -> tuple:
    """Multiprocessing worker task function using the persistent worker handle."""
    global _process_local_al_handle
    return test_damage_with_handle(
        ref=ref,
        al_handle=_process_local_al_handle,
        wlen=wlen,
        g2a=g2a,
        subsample=subsample,
        show_al=show_al,
        verbose=verbose,
        need_read_dict=need_read_dict,
    )


def test_damage(
    ref: str,
    bam: str,
    mode: str,
    wlen: int,
    g2a: bool,
    subsample: float,
    show_al: bool,
    process: int,
    verbose: bool,
) -> tuple:
    """Prepare data and run statistical tests for ancient damage on a reference contig (standalone)."""
    al_handle = pysam.AlignmentFile(bam, mode=mode, threads=1)
    try:
        res = test_damage_with_handle(
            ref=ref,
            al_handle=al_handle,
            wlen=wlen,
            g2a=g2a,
            subsample=subsample,
            show_al=show_al,
            verbose=verbose,
            need_read_dict=True,
        )
        al_handle.close()
        return res
    except Exception:
        al_handle.close()
        return False, {}
