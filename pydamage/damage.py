#!/usr/bin/env python3

import numpy as np
import pysam
from .ctot import ct_al


class ref_alignments(pysam.AlignmentFile):

    def __init__(reference):
        """Constructor of the class

        Args:
            reference (string): a reference from the indexed bam file 
        """
        self.alignments = self.fetch(reference)

    def get_ct(wlen, print_al):
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
                                print_al=print_al)
