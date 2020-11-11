#!/usr/bin/env python3


def damage_al(reference, read_name, ref_name, query, cigartuple, wlen, show_al):
    """Compute CtoT mutations for a single alignment

    Args:
        reference (str): reference sequence
        read_name(str): name of read
        ref_name(str): name of reference
        query (string): query sequence
        cigartuple (tuple): cigar tuple (pysam)
        wlen (int): window length
        print_al (bool): print alignment
    Returns:
        dict : {'CT': [ CtoT pos], 
                'GA': [GtoA pos], 
                'A':[A pos],
                'T':[T pos],
                'G':[G pos],
                'C':[C pos],
                'N':[N pos],
                ' ':[no coverage pos],
                'all': [all pos]}
    """
    r_pos = 0
    q_pos = 0
    r_string = ""
    q_string = ""
    res = ""
    base_trans_counts = {'A': [],
                         'T': [],
                         'G': [],
                         'C': [],
                         'N': [],
                         ' ': [],
                         'all': [],
                         'CT': [],
                         'CC': [],
                         'GA': []}
    for c in cigartuple:
        # [M, =, X] - alignment match (can be a sequence match or mismatch)
        if c[0] in [0, 7, 8]:
            r_string += reference[r_pos:r_pos+c[1]]
            q_string += query[q_pos:q_pos+c[1]]
            r_pos += c[1]
            q_pos += c[1]
        # I, S - Insertion of soft clipping
        if c[0] in [1, 4]:
            r_string += " "*c[1]
            q_string += query[q_pos:q_pos+c[1]]
            r_pos = r_pos
            q_pos += c[1]
        # D - Deletion in query: add space in query string
        if c[0] in [2, 3]:
            r_string += reference[r_pos:r_pos+c[1]]
            q_string += " "*c[1]
            r_pos += c[1]
            q_pos = q_pos

    for i in range(len(r_string)):
        r_char = r_string[i].upper()
        q_char = q_string[i].upper()
        base_trans_counts[q_char].append(i)
        if q_char in ['A', 'T', 'G', 'C']:
            base_trans_counts['all'].append(i)
        if r_char == "C" and q_char == "C":
            base_trans_counts['CC'].append(i)
        if r_char != q_char:
            if r_char == "C" and q_char == "T":
                base_trans_counts['CT'].append(i)
            if r_char == "G" and q_char == "A":
                base_trans_counts['GA'].append(i)

    if show_al:
        res += "R " + r_string + "\t" + ref_name + "\n  "
        for i in range(0, min(len(r_string), len(q_string))):
            if r_string[i] == q_string[i]:
                res += ("|")
            else:
                res += (" ")
        res += "\nQ " + q_string + "\t" + read_name + "\n "
        print(res)

    return(base_trans_counts)
