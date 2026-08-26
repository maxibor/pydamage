#!/usr/bin/env python3

import logging
from os import path, makedirs
import matplotlib.pyplot as plt
from matplotlib import use
import numpy as np
import pandas as pd
from pysam import FastxFile
from statsmodels.nonparametric.smoothers_lowess import lowess

from pydamage.models import damage_model, null_model

plt.set_loglevel("WARNING")
use("Agg")


def bin_plot(csv: str, fasta: str, outdir: str, write_fig: bool = True, **kwargs) -> tuple:
    """Plot ancient DNA damage smile plot for a metagenomic bin.

    Args:
        csv (str): Path to PyDamage results CSV file
        fasta (str): Path to bin FASTA file
        outdir (str): Output directory
        write_fig (bool): If True, save plot PNG to outdir

    Returns:
        tuple: (ct_mean, ct_std, ga_mean, ga_std)
    """
    df = pd.read_csv(csv)
    pydam_refs = set(df["reference"])
    binname = path.splitext(path.basename(fasta))[0]
    contigs = set()

    with FastxFile(fasta) as fh:
        for entry in fh:
            contigs.add(entry.name)

    intersect = contigs.intersection(pydam_refs)
    union = contigs.union(pydam_refs)
    if union != intersect:
        logging.warning(f"{len(union - intersect)} contigs are not shared between fasta and pydamage results")

    df = df[df["reference"].isin(intersect)]
    nb_reads = int(df["nb_reads_aligned"].sum())
    nb_contigs = int(df["reference"].nunique())
    reflen = int(df["reflen"].sum())
    coverage = str(round((df["coverage"] * df["reflen"]).sum() / reflen, 2)) if reflen > 0 else "0"

    df_ct = df.filter(regex="CtoT.*")
    df_ct.columns = df_ct.columns.str.replace("CtoT-", "")
    ct_mean = df_ct.mean(axis=0)
    ct_std = df_ct.std(axis=0)

    df_ga = df.filter(regex="GtoA.*")
    df_ga.columns = df_ga.columns.str.replace("GtoA-", "")
    ga_mean = df_ga.mean(axis=0)
    ga_std = df_ga.std(axis=0)

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10, 5))
    title = f"{fasta}\nNumber of reads: {nb_reads} | Number of contigs: {nb_contigs} | Bin length: {reflen} | Coverage: {coverage}\n"
    fig.suptitle(title, fontsize=11)

    ax1.plot(ct_mean, color="#bd0d45", label="mean")
    ax1.fill_between(
        ct_mean.index, ct_mean - ct_std, ct_mean + ct_std, alpha=0.2, color="#bd0d45", label="± 1 SD"
    )
    ax1.set_xlabel("Position in read from 5'", fontsize=10)
    ax1.set_xticks(ax1.get_xticks())
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45)
    ax1.legend()
    ax1.set_title("C to T substitution frequency", fontsize=8)
    ax1.set_ylabel("Substitution frequency", fontsize=10)

    ax2.set_xlabel("Position in read from 3'", fontsize=10)
    ax2.set_title("G to A substitution frequency", fontsize=8)

    if len(ga_mean) > 0:
        ax2.plot(ga_mean, color="#0000FF", label="mean")
        ax2.fill_between(
            ga_mean.index, ga_mean - ga_std, ga_mean + ga_std, alpha=0.2, color="#0000FF", label="± 1 SD"
        )
        ax2.invert_xaxis()
        ax2.set_xticks(ax2.get_xticks())
        ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45)
        ax2.legend()
    else:
        ax2.set_xticks(ax1.get_xticks())
        ax2.set_xticklabels([None] * len(ax2.get_xticklabels()))

    if write_fig:
        if outdir and not path.exists(outdir):
            makedirs(outdir, exist_ok=True)
        outpath = path.join(outdir, f"{binname}.png") if outdir else f"{binname}.png"
        fig.savefig(outpath, dpi=300)
        logging.info(f"Wrote PyDamage bin plot to {outpath}")

    plt.close(fig)
    return ct_mean, ct_std, ga_mean, ga_std


def damageplot(damage_dict: dict, wlen: int, plot_g2a: bool, outdir: str):
    """Draw and save damage plot with model fit and residuals for a contig."""
    x = np.arange(wlen)
    c2t = np.array([damage_dict.get(f"CtoT-{i}", 0.0) for i in x])
    g2a = np.array([damage_dict.get(f"GtoA-{i}", 0.0) for i in x]) if plot_g2a else None

    p0 = damage_dict["p0"]
    p0_stdev = damage_dict["p0_stdev"]
    p = damage_dict["p"]
    pmin = damage_dict["pmin"]
    pmin_stdev = damage_dict["pmin_stdev"]
    pmax = damage_dict["pmax"]
    pmax_stdev = damage_dict["pmax_stdev"]
    contig = damage_dict["reference"]
    pvalue = damage_dict["pvalue"]
    coverage = damage_dict["coverage"]
    residuals = damage_dict["residuals"]
    rmse = damage_dict["RMSE"]

    rpval = "<0.001" if pvalue < 0.001 else f"={round(pvalue, 3)}"

    m_null = null_model()
    p0_low = max(m_null.bounds[0][0], p0 - 2 * p0_stdev)
    p0_high = min(m_null.bounds[1][0], p0 + 2 * p0_stdev)
    y_unif = m_null.fit(x, p0)
    y_unif_low = np.maximum(np.zeros(y_unif.shape[0]), m_null.fit(x, p0_low))
    y_unif_high = np.minimum(np.ones(y_unif.shape[0]), m_null.fit(x, p0_high))

    geom = damage_model()
    geom_pmin_low = max(geom.bounds[0][1], pmin - 2 * pmin_stdev)
    geom_pmin_high = min(geom.bounds[1][1], pmin + 2 * pmin_stdev)
    geom_pmax_low = max(geom.bounds[0][2], pmax - 2 * pmax_stdev)
    geom_pmax_high = min(geom.bounds[1][2], pmax + 2 * pmax_stdev)

    y_geom = geom.fit(x, p, pmin, pmax, wlen=wlen)
    y_geom_low = np.maximum(np.zeros(y_geom.shape[0]), geom.fit(x, p, geom_pmin_low, geom_pmax_low, wlen=wlen))
    y_geom_high = np.minimum(np.ones(y_geom.shape[0]), geom.fit(x, p, geom_pmin_high, geom_pmax_high, wlen=wlen))

    fig, ax = plt.subplots()

    # Null model
    ax.plot(x, y_unif, linewidth=2.5, color="DarkOliveGreen", alpha=0.8, label="Null model")
    ax.fill_between(x, y_unif_low, y_unif_high, color="DarkOliveGreen", alpha=0.1, label="Null Model CI (2 sigma)")

    # Damage model
    ax.plot(x, y_geom, linewidth=2.5, color="#D7880F", alpha=0.8, label="Damage model")
    ax.fill_between(x, y_geom_low, y_geom_high, color="#D7880F", alpha=0.1, label="Damage Model CI (2 sigma)")

    # Observed transitions
    ax.plot(x, c2t, color="#bd0d45", alpha=0.2, label="C to T transitions from 5' end (forward)")
    if plot_g2a and g2a is not None:
        ax.plot(x, g2a, color="#0000FF", alpha=0.2, label="G to A transitions from 3' end (reverse)")

    ax.set_xlabel("Position in read", fontsize=10)
    ax.set_ylabel("Substitution frequency", fontsize=10)
    step = 5 if wlen >= 10 else 1
    ticks = np.arange(x[0], x[-1] + 1, step)
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticks, rotation=45, fontsize=6)
    ax.set_title(f"coverage: {round(coverage, 2)} - pvalue{rpval}", fontsize=8)
    ax.legend(fontsize=8)

    # Inset: Residuals vs Fitted
    left, bottom, width, height = [0.65, 0.4, 0.2, 0.2]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.set_ylim([-max(pmax, 0.01), max(pmax, 0.01)])
    fitted = x
    smoothed = lowess(residuals, fitted)
    ax2.scatter(fitted, residuals, marker=".", color="black")
    ax2.plot(smoothed[:, 0], smoothed[:, 1], color="r")
    ax2.plot([min(fitted), max(fitted)], [0, 0], color="k", linestyle=":", alpha=0.3)
    ax2.set_ylabel("Residuals", fontsize=6)
    ax2.set_xlabel("Fitted Values", fontsize=6)
    ax2.set_title(f"Residuals vs. Fitted\nRMSE={round(rmse, 3)}", fontsize=6)
    ax2.set_xticks(ticks)
    ax2.set_xticklabels(ticks, fontsize=6, rotation=45)
    fig.suptitle(contig, fontsize=12, y=0.95)

    if outdir and not path.exists(outdir):
        makedirs(outdir, exist_ok=True)
    outpath = path.join(outdir, f"{contig}.png") if outdir else f"{contig}.png"
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
