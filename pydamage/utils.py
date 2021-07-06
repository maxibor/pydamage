import numpy as np
import os
import sys
import shutil
from statsmodels.stats.multitest import multipletests
import pandas as pd
from typing import Tuple


def check_extension(filename: str) -> str:
    """Check alignment file format to give correct open mode

    Args:
        filename (str): Path to alignment file

    Returns:
        str: opening mode

    Raises:
        Exception: Extension not supported
    """
    extension = filename.split(".")[-1]
    modes = {"bam": "rb", "sam": "r", "cram": "rc"}
    try:
        return modes[extension]
    except KeyError:
        raise Exception(f"{extension} file extension not supported")


def makedir(dirpath: str, confirm: bool = True, force: bool = False):
    """Make directory if user confirms overwritting

    Args:
        dirpath (str): Path to directory
        confirm (bool, optional): Ask the user to confirm. Defaults to True.
        force (bool, optional): Always create directory. Defaults to False.
    """
    if os.path.exists(dirpath):
        if confirm and force is False:
            print(
                f"Result directory, {dirpath}, already exists, it will be overwritten"
            )
            if input("Do You Want To Continue? (y|n) ").lower() != "y":
                sys.exit()
        shutil.rmtree(dirpath)

    os.makedirs(dirpath)


def pandas_processing(res_dict: dict, wlen: int) -> pd.core.frame.DataFrame:
    """Performs Pandas processing of Pydamage results

    Args:
        res_dict (dict): Result dictionary of LR test
        wlen(int): window size
    """
    df = pd.DataFrame(res_dict)
    if len(res_dict) == 0:
        return df
    qvalues = pd.Series(
        multipletests(df["pvalue"].dropna(), method="fdr_bh")[1],
        index=df["pvalue"].dropna().index,
        name="qvalue",
    )
    df = df.merge(qvalues, left_index=True, right_index=True, how="outer")
    df = df[
        [
            "p0",
            "p0_stdev",
            "p",
            "p_stdev",
            "pmin",
            "pmin_stdev",
            "pmax",
            "pmax_stdev",
            "pvalue",
            "qvalue",
            "RMSE",
            "reference",
            "nb_reads_aligned",
            "coverage",
            "reflen",
        ]
        + [f"CtoT-{i}" for i in range(wlen)]
    ]
    df.rename(
        columns={
            "p0": "null_model_p0",
            "p0_stdev": "null_model_p0_stdev",
            "p": "damage_model_p",
            "p_stdev": "damage_model_p_stdev",
            "pmin": "damage_model_pmin",
            "pmin_stdev": "damage_model_pmin_stdev",
            "pmax": "damage_model_pmax",
            "pmax_stdev": "damage_model_pmax_stdev",
        },
        inplace=True,
    )
    df.sort_values(by=["qvalue"], inplace=True)
    df.set_index("reference", inplace=True)
    df.dropna(axis=1, how="all", inplace=True)
    df = df.round(3)
    return df


def pandas_group_processing(res_dict: dict) -> pd.core.frame.DataFrame:
    """Performs Pandas processing of Pydamage grouped reference results

    Args:
        res_dict (dict): Result dictionary of LR test
    """
    filt_dict = {}
    filt_dict["reference"] = str(res_dict["reference"])
    filt_dict["null_model_p0"] = float(res_dict["p0"])
    filt_dict["null_model_p0_stdev"] = float(res_dict["p0_stdev"])
    filt_dict["damage_model_p"] = float(res_dict["p"])
    filt_dict["damage_model_p_stdev"] = float(res_dict["p_stdev"])
    filt_dict["damage_model_pmin"] = float(res_dict["pmin"])
    filt_dict["damage_model_pmin_stdev"] = float(res_dict["pmin_stdev"])
    filt_dict["damage_model_pmax"] = float(res_dict["pmax"])
    filt_dict["damage_model_pmax_stdev"] = float(res_dict["pmax_stdev"])
    filt_dict["pvalue"] = float(res_dict["pvalue"])
    filt_dict["RMSE"] = float(res_dict["RMSE"])
    filt_dict["nb_reads_aligned"] = int(res_dict["nb_reads_aligned"])
    filt_dict["coverage"] = float(res_dict["coverage"])
    filt_dict["reflen"] = int(res_dict["reflen"])

    for i in list(res_dict.keys()):
        if str(i).startswith("CtoT"):
            filt_dict[i] = float(res_dict[i])
        if str(i).startswith("GtoA"):
            filt_dict[i] = float(res_dict[i])

    df = (
        pd.Series(filt_dict)
        .to_frame(name="reference")
        .transpose()
        .set_index("reference")
    )
    return df


def df_to_csv(
    df: pd.core.frame.DataFrame, outdir: str, outfile: str = "pydamage_results.csv"
):
    """Write Pydamage results to disk

    Args:
        df(pandas DataFrame): Pydamage results DataFrame
        outdir (str): Path to output directory
    """
    df = df.round(3)
    if not outdir:
        outdir = "."
    df.to_csv(f"{outdir}/{outfile}")


def sort_dict_by_keys(adict: dict) -> dict:
    """Sort dictonary by keys

    Args:
        adict (dict): dictorary to sort
    Returns:
        dict: dictionary sorted on keys
    """
    res = {}
    for k in sorted(adict.keys()):
        res[k] = adict[k]
    return res


def RMSE(residuals: np.ndarray) -> float:
    """Computes Root Mean Square Error

    Args:
        residuals (np.array(float)): Array of residuals
    Returns:
        float: RMSE
    """
    return np.sqrt(np.mean(residuals ** 2))


def create_damage_dict(
    damage_data: list, non_damage_data: list, wlen: int
) -> Tuple[dict, dict]:
    """Creates C bases positions dictionnary

    For C->T transitions
    For C in reference and in query

    Args:
        damage_data (list of int): List of positions where CtoT transitions were observed
        non_damage_data (list of int): List of positions where C in ref and query
    Returns:
        dict{int:int}: {position: number of CtoT transitions at position}
        dict{int:int}: {position: number of C match at position}
    """
    damage_pos, damage_counts = np.unique(np.sort(damage_data), return_counts=True)
    non_damage_pos, non_damage_counts = np.unique(
        np.sort(non_damage_data), return_counts=True
    )
    damage = dict(zip(damage_pos, damage_counts))
    non_damage = dict(zip(non_damage_pos, non_damage_counts))
    damage_dict = {}
    non_damage_dict = {}

    for i in range(wlen):
        if i in non_damage or i in damage:
            if i in damage:
                damage_dict[i] = damage[i]
            else:
                damage_dict[i] = 0
            if i in non_damage:
                non_damage_dict[i] = non_damage[i]
            else:
                non_damage_dict[i] = 0

    return (damage_dict, non_damage_dict)
