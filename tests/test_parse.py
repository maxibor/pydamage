import os
import pytest
from pydamage import parse_ct


@pytest.fixture(autouse=True)
def bam():
    os
    bam = os.path.dirname(os.path.realpath(__file__))+"/data/aligned.bam"
    return(bam)


def test_parsing():
    pass
