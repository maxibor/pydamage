from matplotlib import use
import matplotlib.pyplot as plt
from pydamage.models import damage_model, null_model
from statsmodels.nonparametric.smoothers_lowess import lowess
from pydamage import utils
import numpy as np
from os import makedirs
from scipy.stats import probplot

use("Agg")


def damageplot(damage_dict, outdir):
    """Draw pydamage plots

    Args:
        damage_dict(dict): pydamage result dictionary
        qlen(int): query length
        outdir(str): Pydamage result directory
    """
    x = np.array(range(damage_dict['wlen']))
    qlen = np.array(range(damage_dict['qlen']))
    y = np.array([damage_dict[i] for i in x])
    c2t = np.array([damage_dict[f"CtoT-{i}"] for i in qlen])
    g2a = np.array([damage_dict[f"GtoA-{i}"] for i in qlen])
    p0 = damage_dict['p0']
    p0_stdev = damage_dict['p0_stdev']
    p = damage_dict['p']
    pmin = damage_dict['pmin']
    pmin_stdev = damage_dict['pmin_stdev']
    pmax = damage_dict['pmax']
    pmax_stdev = damage_dict['pmax_stdev']
    contig = damage_dict['reference']
    pvalue = damage_dict['pvalue']
    coverage = damage_dict['coverage']
    residuals = damage_dict['residuals']
    rmse = damage_dict['RMSE']
    plotdir = outdir

    if pvalue < 0.001:
        rpval = "<0.001"
    else:
        rpval = f"={round(pvalue,3)}"

    m_null = null_model()
    p0_low = max(m_null.bounds[0][0], p0 - 2*p0_stdev)
    p0_high = min(m_null.bounds[1][0], p0 + 2*p0_stdev)
    y_unif = m_null.fit(x, p0)
    y_unif_low = np.maximum(
        np.zeros(y_unif.shape[0]), m_null.fit(x, p0_low))
    y_unif_high = np.minimum(
        np.ones(y_unif.shape[0]), m_null.fit(x, p0_high))

    geom = damage_model()
    geom_pmin_low = max(geom.bounds[0][1], pmin - 2*pmin_stdev)
    geom_pmin_high = min(geom.bounds[1][1], pmin + 2*pmin_stdev)
    geom_pmax_low = max(geom.bounds[0][2], pmax - 2*pmax_stdev)
    geom_pmax_high = min(geom.bounds[1][2], pmax + 2*pmax_stdev)

    y_geom = geom.fit(x, p, pmin, pmax)
    y_geom_low = np.maximum(np.zeros(y_geom.shape[0]), geom.fit(
        x, p, geom_pmin_low, geom_pmax_low))
    y_geom_high = np.minimum(np.ones(y_geom.shape[0]), geom.fit(
        x, p, geom_pmin_high, geom_pmax_high))

    plt.xticks(rotation=45, fontsize=8)

    fig, ax = plt.subplots()

    ax.plot(x, y_unif,
            linewidth=2.5,
            color='DarkOliveGreen',
            alpha=0.8,
            label='Null model')

    ax.fill_between(x, y_unif_low, y_unif_high,
                    color='DarkOliveGreen',
                    alpha=0.1,
                    label='Null Model CI (2 sigma)')

    ax.plot(x, y_geom,
            linewidth=2.5,
            color='#D7880F',
            alpha=0.8,
            label='Damage model')

    ax.fill_between(x, y_geom_low, y_geom_high,
                    color='#D7880F',
                    alpha=0.1,
                    label='Damage Model CI (2 sigma)')

    ax.plot(qlen, g2a,
            color='#236cf5',
            alpha=0.1,
            label='G to A transitions')

    ax.plot(qlen, c2t,
            color='#bd0d45',
            alpha=0.2,
            label='C to T transitions')

    ax.set_xlabel("Base from 5'", fontsize=10)
    ax.set_ylabel("Substitution frequency", fontsize=10)
    ax.xaxis.set_ticks(np.arange(qlen[0], qlen[-1], 5))
    ax.set_xticklabels(ax.get_xticks(), rotation=45, fontsize=6)
    ax.set_title(f"coverage: {round(coverage,2)} - pvalue{rpval}", fontsize=8)
    ax.legend(fontsize=8)
    # ax.set_title(f"coverage: {round(coverage,2)} | pvalue{rpval}", fontsize=8)

    left, bottom, width, height = [0.65, 0.3, 0.2, 0.2]
    ax2 = fig.add_axes([left, bottom, width, height])

    fitted = x
    smoothed = lowess(residuals, fitted)
    ax2.scatter(fitted, residuals, marker=".", color='black')
    ax2.plot(smoothed[:, 0], smoothed[:, 1], color='r')
    ax2.plot([min(fitted), max(fitted)], [0, 0],
             color='k', linestyle=':', alpha=.3)
    ax2.set_ylabel('Residuals', fontsize=6)
    ax2.set_xlabel('Fitted Values', fontsize=6)
    ax2.set_title(f"Residuals vs. Fitted\nRMSE={round(rmse, 3)}", fontsize=6)
    ax2.xaxis.set_ticks(np.arange(fitted[0], fitted[-1], 5))
    ax2.set_xticklabels([int(i) for i in ax2.get_xticks()],
                        fontsize=6, rotation=45)
    ax2.set_yticklabels([round(i, 3) for i in ax2.get_yticks()], fontsize=6)

    fig.suptitle(contig, fontsize=12, y=0.95)

    fig.savefig(f"{plotdir}/{contig}.png", dpi=200)
