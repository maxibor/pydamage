import pytest
import pysam
import random
import numpy as np


@pytest.fixture(autouse=True)
def bamfile():
    return pysam.AlignmentFile("tests/data/aligned.bam", "rb")


def test_al_to_damage(bamfile):
    from pydamage.damage import al_to_damage

    al = al_to_damage(reference="NZ_JHCB02000002.1", al_handle=bamfile)
    all_ct, all_ga, all_cc, all_c, all_g, all_bases = al.get_damage(
        wlen=20, show_al=False
    )

    assert all_ct == [
        15,
        46,
        74,
        0,
        0,
        42,
        2,
        82,
        66,
        15,
        43,
        22,
        11,
        39,
        55,
        0,
        37,
        38,
    ]
    assert all_ga == [
        55,
        1,
        64,
        42,
        0,
        27,
        6,
        41,
        25,
        28,
        65,
        40,
        102,
        49,
        68,
        19,
        10,
        43,
        13,
    ]
    assert all_cc[:10] == [7, 11, 18, 39, 43, 45, 50, 51, 52, 63]
    assert all_c[:10] == [7, 11, 18, 39, 43, 45, 50, 51, 52, 63]
    assert all_g[:10] == [0, 2, 13, 17, 24, 27, 32, 33, 46, 49]
    random.seed(42)
    assert random.sample(all_bases, 10) == [33, 23, 6, 63, 28, 21, 3, 39, 67, 18]


def test_avg_coverage():
    from pydamage.damage import avg_coverage

    cov = np.array([[1, 4, 5, 3], [2, 5, 6, 7], [5, 6, 3, 1], [2, 4, 8, 1]], np.int32)
    assert avg_coverage(cov) == 15.75


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
    test_dict3 = {
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
        "base_cov": np.array([153, 153, 153, 0, 153, 153, 153, 153, 153, 153]),
    }

    assert check_model_fit(test_dict, wlen=11, verbose=False) is False
    assert check_model_fit(test_dict, wlen=10, verbose=False) == test_dict
    assert check_model_fit(test_dict2, wlen=10, verbose=False) is False
    assert check_model_fit(test_dict3, wlen=10, verbose=False) is False


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
    assert dam[0] == pytest.approx(0.16666666666666666)
    assert dam["CtoT-0"] == pytest.approx(0.16666666666666666)

    assert dam["p"] == pytest.approx(0.9899999999999994)
    assert dam["pmin"] == pytest.approx(0.011620247207355617)
    assert dam["pmax"] == pytest.approx(0.16653935252468358)
    assert dam["p_stdev"] == pytest.approx(0.20326991191464713)
    assert dam["pmin_stdev"] == pytest.approx(0.00722931975363037)
    assert dam["pmax_stdev"] == pytest.approx(0.030636963039075888)
    assert dam["p0"] == pytest.approx(0.019444463451190847)
    assert dam["p0_stdev"] == pytest.approx(0.010096191035450966)
    assert dam["pvalue"] == pytest.approx(0.01591779796235382)

    assert (
        dam["base_cov"][:10]
        == np.array([153, 153, 153, 153, 153, 153, 153, 153, 153, 153])
    ).all()
    assert dam["model_params"] == pytest.approx(
        [
            0.9899999999999994,
            0.011620247207355617,
            0.16653935252468358,
            0.019444463451190847,
            0.20326991191464713,
            0.00722931975363037,
            0.030636963039075888,
            0.010096191035450966,
        ]
    )
    assert dam["qlen"] == 130
    assert dam["RMSE"] == pytest.approx(0.028245884410886473)
    assert dam["wlen"] == 20
    assert dam["reference"] == "NZ_JHCB02000002.1"
    assert dam["nb_reads_aligned"] == 153
    assert dam["coverage"] == pytest.approx(0.18994194094919317)
    assert dam["reflen"] == 48399


@pytest.fixture(autouse=True)
def test_get_damage_group():
    from pydamage.damage import get_damage_group

    dmg_grp = get_damage_group(
        ref="NZ_JHCB02000002.1",
        bam="tests/data/aligned.bam",
        mode="rb",
        show_al=False,
        wlen=20,
        process=1,
    )

    # ct_data
    assert dmg_grp[0] == [
        15,
        46,
        74,
        0,
        0,
        42,
        2,
        82,
        66,
        15,
        43,
        22,
        11,
        39,
        55,
        0,
        37,
        38,
    ]
    # ga_data
    assert dmg_grp[1] == [
        55,
        1,
        64,
        42,
        0,
        27,
        6,
        41,
        25,
        28,
        65,
        40,
        102,
        49,
        68,
        19,
        10,
        43,
        13,
    ]
    # cc_data
    assert dmg_grp[2][:10] == [7, 11, 18, 39, 43, 45, 50, 51, 52, 63]
    #  c_data
    assert dmg_grp[3][:10] == [7, 11, 18, 39, 43, 45, 50, 51, 52, 63]
    # g_data
    assert dmg_grp[4][:10] == [0, 2, 13, 17, 24, 27, 32, 33, 46, 49]
    # all_bases
    random.seed(42)
    assert random.sample(dmg_grp[5], 10) == [33, 23, 6, 63, 28, 21, 3, 39, 67, 18]
    # cov
    assert dmg_grp[6] == pytest.approx(0.18994194094919317)
    # nb_reads_aligned
    assert dmg_grp[7] == 153
    # reflen
    assert dmg_grp[8] == 48399

    return dmg_grp


def test_test_damage_group(test_get_damage_group):
    from pydamage.damage import test_damage_group

    damg = test_damage_group(
        ct_data=test_get_damage_group[0],
        ga_data=test_get_damage_group[1],
        cc_data=test_get_damage_group[2],
        all_bases=test_get_damage_group[5],
        nb_reads_aligned=test_get_damage_group[7],
        cov=test_get_damage_group[6],
        reflen=test_get_damage_group[8],
        wlen=20,
        verbose=False,
    )

    assert damg[0] == pytest.approx(0.16666666666666666, rel=1e-3)
    assert damg[10] == pytest.approx(0.0)
    assert damg["p"] == pytest.approx(0.9899999999999994, rel=1e-3)
    assert damg["pmin"] == pytest.approx(0.011620247207355617, rel=1e-3)
    assert damg["pmax"] == pytest.approx(0.16653935252468358, rel=1e-3)
    assert damg["p_stdev"] == pytest.approx(0.20326991191464713, rel=1e-3)
    assert damg["pmin_stdev"] == pytest.approx(0.00722931975363037, rel=1e-3)
    assert damg["pmax_stdev"] == pytest.approx(0.030636963039075888, rel=1e-3)
    assert damg["p0"] == pytest.approx(0.019444463451190847, rel=1e-3)
    assert damg["p0_stdev"] == pytest.approx(0.010096191035450966, rel=1e-3)
    assert damg["pvalue"] == pytest.approx(0.01591779796235382, rel=1e-3)

    assert (
        damg["base_cov"][:10]
        == np.array([153, 153, 153, 153, 153, 153, 153, 153, 153, 153])
    ).all()

    assert damg["model_params"][:4] == pytest.approx(
        [
            0.9899999999999994,
            0.011620247207355617,
            0.16653935252468358,
            0.019444463451190847,
        ]
    )
    assert damg["RMSE"] == pytest.approx(0.028245884410886473)
    assert damg["wlen"] == 20
    assert damg["reference"] == "reference"
    assert damg["nb_reads_aligned"] == 153
    assert damg["coverage"] == pytest.approx(0.18994194094919317)
    assert damg["reflen"] == 48399
