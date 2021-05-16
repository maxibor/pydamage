#!/usr/bin/env python3


def damage_al(
    reference,
    read_name,
    is_reverse,
    ref_name,
    query,
    cigartuple,
    wlen,
    show_al,
):
    """Compute CtoT mutations for a single alignment

    Args:
        reference (str): reference sequence
        read_name(str): name of read
        is_reverse(bool): alignment if reverse
        ref_name(str): name of reference
        query (string): query sequence
        cigartuple (tuple): cigar tuple (pysam)
        wlen (int): window length
        print_al (bool): print alignment
    Returns:
        dict : {'C':  [ C pos from 5'],
                'CT': [ CtoT pos from 5'],
                'G':  [ G pos from 3'],
                'GA': [GtoA pos from 3'],
                'no_mut': [Positions where C , and G (if reversed) are conserved]}
    """
    r_pos = 0
    q_pos = 0
    r_string = ""
    q_string = ""
    res = ""
    base_trans_counts = {"C": [], "CT": [], "no_mut": []}
    for c in cigartuple:
        # [M, =, X] - alignment match (can be a sequence match or mismatch)
        if c[0] in [0, 7, 8]:
            r_string += reference[r_pos : r_pos + c[1]]
            q_string += query[q_pos : q_pos + c[1]]
            r_pos += c[1]
            q_pos += c[1]
        # I, S - Insertion of soft clipping
        if c[0] in [1, 4]:
            r_string += " " * c[1]
            q_string += query[q_pos : q_pos + c[1]]
            r_pos = r_pos
            q_pos += c[1]
        # D - Deletion in query: add space in query string
        if c[0] in [2, 3]:
            r_string += reference[r_pos : r_pos + c[1]]
            q_string += " " * c[1]
            r_pos += c[1]
            q_pos = q_pos

    read_len = len(r_string)
    for i in range(min(read_len, wlen)):
        r_char = r_string[i].upper()
        q_char = q_string[i].upper()
        # base_trans_counts[q_char].append(i)
        if is_reverse is False:
            if r_char == "C":
                base_trans_counts["C"].append(i)
                if q_char == "T":
                    base_trans_counts["CT"].append(i)
                elif q_char == "C":
                    base_trans_counts["no_mut"].append(i)
    if show_al:
        orient = {False: ["5'", "3'"], True: ["3'", "5'"]}
        res += "R   " + r_string + "\t" + ref_name + "\n" + orient[is_reverse][0] + "  "
        for i in range(0, min(len(r_string), len(q_string))):
            if r_string[i] == q_string[i]:
                res += "|"
            else:
                res += " "
        res += orient[is_reverse][1] + "\nQ   " + q_string + "\t" + read_name + "\n "
        print(res)

    return base_trans_counts
