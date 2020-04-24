#!/usr/bin/env python3


def damage_al(reference, query, cigartuple, wlen, show_al):
    """Compute CtoT mutations for a single alignment

    Args:
        reference (string): reference sequence
        query (string): query sequence
        cigartuple (tuple): cigar tuple (pysam)
        wlen (int): window length
        print_al (bool): print alignment
    Returns:
        CT (list): list of CtoT positions
        GA (list): list of GtoA positions
    """
    r_pos = 0
    q_pos = 0
    r_string = ""
    q_string = ""
    res = ""
    CT = []
    GA = []
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

    for i in range(0, min(wlen, len(r_string))):
        r_char = r_string[i].upper()
        q_char = q_string[i].upper()
        if r_char != q_char:
            if r_char == "C" and q_char == "T":
                CT.append(i)
            if r_char == "G" and q_char == "A":
                GA.append(i)

    if show_al:
        res += "R " + r_string + "\n  "
        for i in range(0, min(len(r_string), len(q_string))):
            if r_string[i] == q_string[i]:
                res += ("|")
            else:
                res += (" ")
        res += "\nQ " + q_string
        print(res)

    return({'CT':CT,'GA':GA})
