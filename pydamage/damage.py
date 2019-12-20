#!/usr/bin/env python3

import numpy as np
import pysam
from parse_ct import ct_al
from vuong import vuong_closeness
import utils
import models


class al_to_ct():

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

    def get_ct(self, wlen, show_al):
        """Compute CtoT substitutions

        Args:
            wlen (int): window length
            print_al(bool): print alignments representations
        """

        all_ct = []
        for al in self.alignments:
            if al.is_reverse == False and al.is_unmapped == False:
                cigar = al.cigartuples
                ref = al.get_reference_sequence()
                quer = al.query_sequence

                all_ct += ct_al(reference=ref,
                                query=quer,
                                cigartuple=cigar,
                                wlen=wlen,
                                show_al=show_al)
        return(all_ct)


def test_ct(ref, bam, mode, wlen, show_al, min_al, process):
    al_handle = pysam.AlignmentFile(bam, mode=mode, threads=process)
    if al_handle.count(contig=ref) > min_al:
        al = al_to_ct(reference=ref, al_handle=al_handle)
        ct_data = al.get_ct(wlen=wlen, show_al=show_al)
        if ct_data:
            model_A = models.unif_mod()
            model_B = models.geom_mod()
            print(f"=====\n{ref}\n")
            test_res = vuong_closeness(
                model_A=model_A, model_B=model_B, data=ct_data)
            print(f"\n=====")
            test_res['reference'] = ref
            return(test_res)
    else:
        pass
