import matplotlib.pyplot as plt
from pydamage.models import geom_mod, unif_mod
import numpy as np

def extract_plot_data(filt_res_dict, wlen):
    res = {}
    res['x'] = np.array(range(wlen))
    res['y'] = np.array([filt_res_dict[i] for i in res['x']])
    res['unif_pmin'] = filt_res_dict['unif_pmin']
    res['unif_pmin_stdev'] = filt_res_dict['unif_pmin_stdev']
    res['geom_p'] = filt_res_dict['geom_p']
    res['geom_p_stdev'] = filt_res_dict['geom_p_stdev']
    res['geom_pmin'] = filt_res_dict['geom_pmin']
    res['geom_pmin_stdev'] = filt_res_dict['geom_pmin_stdev']
    res['geom_pmax'] = filt_res_dict['geom_pmax']
    res['geom_pmax_stdev'] = filt_res_dict['geom_pmax_stdev']
    res['contig'] = filt_res_dict['reference']
    return(res)


def plot_damage(x, 
                y, 
                unif_pmin, 
                unif_pmin_stdev, 
                geom_p, 
                geom_p_stdev,
                geom_pmin, 
                geom_pmin_stdev,
                geom_pmax, 
                geom_pmax_stdev,
                contig,
                outdir='.',
                **kwargs):

    unif = unif_mod()
    unif_pmin_low = max(unif.bounds[0][0], unif_pmin - 2*unif_pmin_stdev)
    unif_pmin_high = min(unif.bounds[1][0], unif_pmin + 2*unif_pmin_stdev)
    y_unif = unif.pmf(x, unif_pmin)
    y_unif_low = np.maximum(np.zeros(y_unif.shape[0]), unif.pmf(x, unif_pmin_low)) 
    y_unif_high = np.minimum(np.ones(y_unif.shape[0]), unif.pmf(x, unif_pmin_high))

    geom = geom_mod()
    # geom_p_low = max(geom.bounds[0][0], geom_p - 2*geom_p_stdev)
    # geom_p_high = min(geom.bounds[1][0], geom_p + 2*geom_p_stdev)
    # geom_pmin_low = max(geom.bounds[0][1], geom_pmin - 2*geom_pmin_stdev)
    # geom_pmin_high = min(geom.bounds[1][1], geom_pmin + 2*geom_pmin_stdev)
    # geom_pmax_low = max(geom.bounds[0][2], geom_pmax - 2*geom_pmax_stdev)
    # geom_pmax_high = min(geom.bounds[1][2], geom_pmax + 2*geom_pmax_stdev)

    geom_p_low = geom_p - 2*geom_p_stdev
    geom_p_high = geom_p + 2*geom_p_stdev
    geom_pmin_low = geom_pmin - 2*geom_pmin_stdev
    geom_pmin_high = geom_pmin + 2*geom_pmin_stdev
    geom_pmax_low = geom_pmax - 2*geom_pmax_stdev
    geom_pmax_high = geom_pmax + 2*geom_pmax_stdev

    y_geom = geom.pmf(x, geom_p, geom_pmin, geom_pmax)
    # y_geom_low = np.maximum(np.zeros(y_geom.shape[0]), geom.pmf(x, geom_p_low, geom_pmin_low, geom_pmax_low))
    # y_geom_high = np.minimum(np.ones(y_geom.shape[0]), geom.pmf(x, geom_p_high, geom_pmin_high, geom_pmax_high))
    y_geom_low = geom.pmf(x, geom_p_low, geom_pmin_low, geom_pmax_low)
    y_geom_high = geom.pmf(x, geom_p_high, geom_pmin_high, geom_pmax_high)


    fig = plt.figure(figsize=(12,8), dpi=100, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    ax.xaxis.labelpad = 20
    ax.yaxis.labelpad = 20

    plt.plot(x, y,
             'o',
             label='Observed damage')
    # plt.hold(True)

    plt.plot(x, y_unif, 
        linewidth=2.5, 
        color = 'IndianRed',
        alpha = 0.8,
        label = 'Uniform model'
    )
    plt.fill_between(x, y_unif_low, y_unif_high,
        color='IndianRed',
        alpha=0.1,
        label = 'Uniform CI (2 sigma)'
    )

    plt.plot(x, y_geom,
             linewidth=2.5, 
             color = 'DarkOliveGreen',
             alpha = 0.8,
             label = 'Geometric model'
    )
    plt.fill_between(x, y_geom_low, y_geom_high,
        color='DarkOliveGreen',
        alpha=0.1,
        label = 'Geometric CI (2 sigma)'
    )
    plt.xlabel("Position", fontsize=20)
    plt.ylabel("Damage", fontsize=20)
    plt.title(contig, fontsize=20)
    ax.legend(fontsize=18)
    ax.set_xticks(x)
    plt.savefig(f"{outdir}/{contig}.png", dpi=200)