import pytest
import numpy as np
import pysam
from pydamage.parse_damage import damage_al
from pydamage.damage import al_to_damage, sort_count_array_dict
from pydamage import utils
from pydamage.rescale import phred_to_prob, compute_new_prob, rescale_qual
from pydamage.exceptions import AlignmentFileError


def test_damage_al_forward_c_to_t():
    # Forward strand alignment with C->T transition at position 0
    ref = "CGTACG"
    query = "TGTACG"
    cigar = [(0, 6)]  # 6M
    res = damage_al(
        reference=ref,
        read_name="read1",
        is_reverse=False,
        ref_name="ref1",
        query=query,
        cigartuple=cigar,
        wlen=6,
        g2a=True,
        show_al=False,
    )
    assert res["C"] == [0, 4]  # positions of C in ref (0 and 4)
    assert res["CT"] == [0]    # C->T mutation at pos 0
    assert res["no_mut"] == [4] # Conserved C at pos 4
    assert res["G"] == []
    assert res["GA"] == []


def test_damage_al_reverse_g_to_a():
    # Reverse strand alignment with G->A transition from 3' end
    ref = "CGTACG"
    query = "CGTACA"  # last base is G->A in reverse context
    cigar = [(0, 6)]
    res = damage_al(
        reference=ref,
        read_name="read2",
        is_reverse=True,
        ref_name="ref1",
        query=query,
        cigartuple=cigar,
        wlen=6,
        g2a=True,
        show_al=False,
    )
    # Ref reversed: G C A T G C (indices 0..5 from 3' end)
    # Query reversed: A C A T G C
    # At pos 0: Ref G, Query A -> G->A mutation at pos 0
    assert 0 in res["GA"]
    assert 0 in res["G"]
    assert 4 in res["G"]
    assert 4 in res["no_mut"]


def test_damage_al_with_indels_and_clipping():
    # Insertion in query (I: 1) and soft clipping (S: 4)
    ref = "ACGTACGT"
    query = "NNACGTACGT"
    cigar = [(4, 2), (0, 8)]  # 2S 8M
    res = damage_al(
        reference=ref,
        read_name="read3",
        is_reverse=False,
        ref_name="ref1",
        query=query,
        cigartuple=cigar,
        wlen=8,
        g2a=False,
        show_al=False,
    )
    assert isinstance(res, dict)
    assert "C" in res
    assert "CT" in res


def test_sort_count_array_dict():
    arr = np.array([3, 1, 3, 2, 1, 3], dtype="uint32")
    pos, counts = sort_count_array_dict(arr)
    np.testing.assert_array_equal(pos, [1, 2, 3])
    np.testing.assert_array_equal(counts, [2, 1, 3])


def test_check_extension():
    assert utils.check_extension("sample.bam") == "rb"
    assert utils.check_extension("sample.sam") == "r"
    assert utils.check_extension("sample.cram") == "rc"
    with pytest.raises(Exception):
        utils.check_extension("sample.invalid")


def test_rescaling_functions():
    # phred_to_prob
    qual = np.array([10, 20, 30], dtype=np.int64)
    probs = phred_to_prob(qual)
    np.testing.assert_allclose(probs, [0.1, 0.01, 0.001], rtol=1e-5)

    # compute_new_prob
    e = np.array([0.01, 0.01])
    d = np.array([0.2, 0.0])
    new_qual = compute_new_prob(e, d)
    assert len(new_qual) == 2
    assert new_qual[0] < 20  # increased error prob reduces Phred score

    # rescale_qual
    read_qual = [30, 30, 30]
    dmg_pmf = np.array([0.2, 0.1, 0.05])
    damage_bases = [0]
    rescaled = rescale_qual(read_qual, dmg_pmf, damage_bases, reverse=False)
    assert len(rescaled) == 3
    assert rescaled[0] < read_qual[0]
