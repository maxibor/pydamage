#!/usr/bin/env python3


def damage_al(
    reference: str,
    read_name: str,
    is_reverse: bool,
    ref_name: str,
    query: str,
    cigartuple: list,
    wlen: int,
    g2a: bool,
    show_al: bool,
) -> dict:
    """Compute CtoT and GtoA mutations for a single alignment.

    Args:
        reference (str): Reference sequence
        read_name (str): Name of read
        is_reverse (bool): True if alignment is reverse complemented
        ref_name (str): Name of reference contig
        query (str): Query sequence
        cigartuple (list): CIGAR tuples (pysam)
        wlen (int): Window length
        g2a (bool): Track G to A transitions
        show_al (bool): Print alignment representation

    Returns:
        dict: {'C': [C positions from 5'],
               'CT': [C->T positions from 5'],
               'G': [G positions from 3'],
               'GA': [G->A positions from 3'],
               'no_mut': [Positions where C, and G (if reverse) are conserved]}
    """
    r_pos = 0
    q_pos = 0
    r_parts = []
    q_parts = []

    for op, length in cigartuple:
        # [M, =, X] - alignment match (sequence match or mismatch)
        if op in (0, 7, 8):
            r_parts.append(reference[r_pos : r_pos + length])
            q_parts.append(query[q_pos : q_pos + length])
            r_pos += length
            q_pos += length
        # I, S - Insertion or soft clipping
        elif op in (1, 4):
            r_parts.append(" " * length)
            q_parts.append(query[q_pos : q_pos + length])
            q_pos += length
        # D, N - Deletion in query or reference skip
        elif op in (2, 3):
            r_parts.append(reference[r_pos : r_pos + length])
            q_parts.append(" " * length)
            r_pos += length

    r_string = "".join(r_parts)
    q_string = "".join(q_parts)
    read_len = len(r_string)
    limit = min(read_len, wlen)

    base_trans_counts = {"C": [], "CT": [], "G": [], "GA": [], "no_mut": []}

    if not is_reverse:
        for i in range(limit):
            r_char = r_string[i].upper()
            q_char = q_string[i].upper()
            if r_char == "C":
                base_trans_counts["C"].append(i)
                if q_char == "T":
                    base_trans_counts["CT"].append(i)
                elif q_char == "C":
                    base_trans_counts["no_mut"].append(i)
    elif g2a:
        for i in range(limit):
            r_char = r_string[read_len - 1 - i].upper()
            q_char = q_string[read_len - 1 - i].upper()
            if r_char == "G":
                base_trans_counts["G"].append(i)
                if q_char == "A":
                    base_trans_counts["GA"].append(i)
                elif q_char == "G":
                    base_trans_counts["no_mut"].append(i)

    if show_al:
        orient = {False: ["5'", "3'"], True: ["3'", "5'"]}
        match_line = "".join(
            "|" if r_string[i] == q_string[i] else " "
            for i in range(min(len(r_string), len(q_string)))
        )
        res = (
            f"R   {r_string}\t{ref_name}\n"
            f"{orient[is_reverse][0]}  {match_line}{orient[is_reverse][1]}\n"
            f"Q   {q_string}\t{read_name}\n "
        )
        print(res)

    return base_trans_counts
