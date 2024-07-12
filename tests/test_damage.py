import pytest
import pysam
import random
import numpy as np


@pytest.fixture(autouse=True)
def bamfile():
    return pysam.AlignmentFile("tests/data/aligned.bam", "rb")


def test_al_to_damage(bamfile):
    from pydamage.damage import al_to_damage

    al = al_to_damage(reference="NZ_JHCB02000002.1", al_handle=bamfile, wlen=20)
    al.get_damage(show_al=False)

    assert al.C[:10] == [7, 11, 18, 2, 18, 9, 7, 11, 15, 4]
    assert al.CT == [15, 0, 0, 2, 15, 11, 0]
    assert al.GA == []
    assert al.damage_bases == [15, 0, 0, 2, 15, 11, 0]
    assert al.C_G_bases[:10] == [7, 11, 18, 2, 18, 9, 7, 11, 15, 4]
    assert al.no_mut[:10] == [7, 11, 18, 2, 18, 9, 7, 11, 4, 16]


def test_avg_coverage():
    from pydamage.damage import avg_coverage_contig

    cov = np.array([[1, 4, 5, 3], [2, 5, 6, 7], [5, 6, 3, 1], [2, 4, 8, 1]], np.int32)
    assert avg_coverage_contig(cov) == 15.75


def test_check_model_fit():
    from pydamage.damage import check_model_fit

    test_dict = {
        "model_params": [
            0.9899999999999994,
            0.011620247207355617,
            0.16653935252468358,
            0.019444463451190847,
            0.20326991191464713,
            0.00722931975363037,
            0.030636963039075888,
            0.010096191035450966,
        ],
        "base_cov": np.array([153, 153, 153, 153, 153, 153, 153, 153, 153, 153]),
    }

    test_dict2 = {
        "model_params": [
            0.9899999999999994,
            0.011620247207355617,
            0.16653935252468358,
            np.inf,
            0.20326991191464713,
            0.00722931975363037,
            0.030636963039075888,
            0.010096191035450966,
        ],
        "base_cov": np.array([153, 153, 153, 153, 153, 153, 153, 153, 153, 153]),
    }

    assert check_model_fit(test_dict, wlen=10, verbose=False) == test_dict
    assert check_model_fit(test_dict2, wlen=10, verbose=False) is False


def test_test_damage():
    from pydamage.damage import test_damage

    dam = test_damage(
        ref="NZ_JHCB02000002.1",
        bam="tests/data/aligned.bam",
        mode="rb",
        show_al=False,
        wlen=20,
        process=1,
        verbose=False,
    )
    assert dam["CtoT-0"] == pytest.approx(0.13043478260869565)

    assert dam["p"] == pytest.approx(0.9899999999803224)
    assert dam["pmin"] == pytest.approx(0.007181450310122387)
    assert dam["pmax"] == pytest.approx(0.1303544631199724)
    assert dam["p_stdev"] == pytest.approx(0.15225306745745634)
    assert dam["pmin_stdev"] == pytest.approx(0.004305275722314164)
    assert dam["pmax_stdev"] == pytest.approx(0.018245226062185038)
    assert dam["p0"] == pytest.approx(0.01340230967468952)
    assert dam["p0_stdev"] == pytest.approx(0.007265315749990326)
    assert dam["pvalue"] == pytest.approx(0.013003463710602237)
    assert dam["model_params"] == pytest.approx(
        [
            0.9899999999803224,
            0.007181450310122387,
            0.1303544631199724,
            0.01340230967468952,
            0.15225306745745634,
            0.004305275722314164,
            0.018245226062185038,
            0.007265315749990326,
        ]
    )
    assert dam["RMSE"] == pytest.approx(0.016821267360790423)
    assert dam["reference"] == "NZ_JHCB02000002.1"
    assert dam["nb_reads_aligned"] == 153
    assert dam["coverage"] == pytest.approx(0.18994194094919317)
    assert dam["reflen"] == 48399
