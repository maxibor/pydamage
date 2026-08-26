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
    show_al: bool = False,
    C_list: list = None,
    CT_list: list = None,
    G_list: list = None,
    GA_list: list = None,
    no_mut_list: list = None,
    ct_ga_out: list = None,
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
        C_list (list): Optional in-place accumulator for C positions
        CT_list (list): Optional in-place accumulator for C->T positions
        G_list (list): Optional in-place accumulator for G positions
        GA_list (list): Optional in-place accumulator for G->A positions
        no_mut_list (list): Optional in-place accumulator for conserved positions
        ct_ga_out (list): Optional list to record positions for read_dict

    Returns:
        dict: {'C': [...], 'CT': [...], 'G': [...], 'GA': [...], 'no_mut': [...]}
    """
    return_dict = C_list is None
    if return_dict:
        C_list = []
        CT_list = []
        G_list = []
        GA_list = []
        no_mut_list = []

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

    if not is_reverse:
        for i in range(limit):
            rc = r_string[i].upper()
            qc = q_string[i].upper()
            if rc == "C":
                C_list.append(i)
                if qc == "T":
                    CT_list.append(i)
                    if ct_ga_out is not None:
                        ct_ga_out.append(i)
                elif qc == "C":
                    no_mut_list.append(i)
    elif g2a:
        for i in range(limit):
            rc = r_string[read_len - 1 - i].upper()
            qc = q_string[read_len - 1 - i].upper()
            if rc == "G":
                G_list.append(i)
                if qc == "A":
                    GA_list.append(i)
                    if ct_ga_out is not None:
                        ct_ga_out.append(i)
                elif qc == "G":
                    no_mut_list.append(i)

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

    if return_dict:
        return {
            "C": C_list,
            "CT": CT_list,
            "G": G_list,
            "GA": GA_list,
            "no_mut": no_mut_list,
        }
    return None
