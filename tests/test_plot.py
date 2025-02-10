import os
import pytest
from pydamage.plot import bin_plot


@pytest.fixture(autouse=True)
def bam():
    os
    bam = os.path.dirname(os.path.realpath(__file__))+"/data/aligned.bam"
    return(bam)


def test_binplot():
    csv = "tests/data/pydamage_results_sequence.csv"
    fasta = "tests/data/sequence.fa"
    ct_mean, ct_std, ga_mean, ga_std = bin_plot(csv, fasta, None, write_fig = False)
    assert ct_mean.to_list()[:3] == pytest.approx([0.033097, 0.032645, 0.015355], abs=1e-3)