#!/usr/bin/env python3

import pysam
from pydamage.parse_damage import damage_al
from pydamage.vuong import vuong_closeness
from pydamage import models
import numpy as np

class al_to_damage():

    def __init__(self, reference, al_handle):
        """Constructor of the class

        Args:
            al_handle(pysam.AlignmentFile)
            mode (str): opening mode (r, rb, rc)
            reference (string): a reference from the indexed bam file 
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
            list: all_ct - positions of reads with CtoT transitions
            list: all_ga - positions of reads with GtoA transitions
        """

        all_ct = []
        all_ga = []
        for al in self.alignments:
            if al.is_reverse == False and al.is_unmapped == False:
                cigar = al.cigartuples
                ref = al.get_reference_sequence()
                quer = al.query_sequence

                all_damage = damage_al(reference=ref,
                                       query=quer,
                                       cigartuple=cigar,
                                       wlen=wlen,
                                       show_al=show_al)
                all_ct += all_damage['CT']
                all_ga += all_damage['GA']

        return(all_ct, all_ga)


def avg_coverage(pysam_cov):
    A = np.array(pysam_cov[0], dtype=int)
    C = np.array(pysam_cov[1], dtype=int)
    G = np.array(pysam_cov[2], dtype=int)
    T = np.array(pysam_cov[3], dtype=int)
    cov_all_bases = A + C + G + T
    cov = np.mean(cov_all_bases)
    return(cov)

def test_damage(ref, bam, mode, wlen, show_al, min_al, min_cov, process, verbose):
    al_handle = pysam.AlignmentFile(bam, mode=mode, threads=process)
    try:
        cov = avg_coverage(al_handle.count_coverage(contig=ref))
        nb_reads_aligned = al_handle.count(contig=ref)
        
        if nb_reads_aligned >= min_al or cov >= min_cov:
            al = al_to_damage(reference=ref, al_handle=al_handle)
            ct_data, ga_data = al.get_damage(wlen=wlen, show_al=show_al)
            if ct_data:
                model_A = models.unif_mod()
                model_B = models.geom_mod()
                test_res = vuong_closeness(ref=ref, 
                                           model_A=model_A, 
                                           model_B=model_B, 
                                           ct_data=ct_data, 
                                           ga_data=ga_data,
                                           wlen=wlen, 
                                           verbose=verbose)
                test_res['reference'] = ref
                test_res['nb_reads_aligned'] = nb_reads_aligned
                test_res['coverage'] = cov
                return(test_res)
        else:
            pass
    except ValueError:
        print(f"Could not fit a model for {ref} because of too few reads aligned ({nb_reads_aligned})")
        pass
