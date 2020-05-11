import matplotlib.pyplot as plt
from pydamage.models import geom_mod, unif_mod
from pydamage import utils
import numpy as np
from os import makedirs


def damageplot(damage_dict, wlen, outdir):
        """Draw pydamage plots
        
        Args:
            damage_dict(dict): pydamage result dictionary
            wlen(int): window length
            qlen(int): query length
            outdir(str): Pydamage result directory
        """        
        x = np.array(range(wlen))
        qlen = np.array(range(damage_dict['qlen']))
        y = np.array([damage_dict[i] for i in x])
        c2t = np.array([damage_dict[f"CtoT-{i}"] for i in qlen])
        g2a = np.array([damage_dict[f"GtoA-{i}"] for i in qlen])
        unif_pmin = damage_dict['unif_pmin']
        unif_pmin_stdev = damage_dict['unif_pmin_stdev']
        geom_p = damage_dict['geom_p']
        geom_pmin = damage_dict['geom_pmin']
        geom_pmin_stdev = damage_dict['geom_pmin_stdev']
        geom_pmax = damage_dict['geom_pmax']
        geom_pmax_stdev = damage_dict['geom_pmax_stdev']
        contig = damage_dict['reference']
        pvalue = damage_dict['pvalue']
        coverage = damage_dict['coverage']
        plotdir = outdir

        if pvalue < 0.001:
            rpval = "<0.001"
        else:
            rpval = f"={round(pvalue,3)}"


        unif = unif_mod()
        unif_pmin_low = max(unif.bounds[0][0], unif_pmin - 2*unif_pmin_stdev)
        unif_pmin_high = min(unif.bounds[1][0], unif_pmin + 2*unif_pmin_stdev)
        y_unif = unif.pmf(x, unif_pmin)
        y_unif_low = np.maximum(np.zeros(y_unif.shape[0]), unif.pmf(x, unif_pmin_low)) 
        y_unif_high = np.minimum(np.ones(y_unif.shape[0]), unif.pmf(x, unif_pmin_high))

        geom = geom_mod()
        geom_pmin_low = max(geom.bounds[0][1], geom_pmin - 2*geom_pmin_stdev)
        geom_pmin_high = min(geom.bounds[1][1], geom_pmin + 2*geom_pmin_stdev)
        geom_pmax_low = max(geom.bounds[0][2], geom_pmax - 2*geom_pmax_stdev)
        geom_pmax_high = min(geom.bounds[1][2], geom_pmax + 2*geom_pmax_stdev)

        y_geom = geom.pmf(x, geom_p, geom_pmin, geom_pmax)
        y_geom_low = np.maximum(np.zeros(y_geom.shape[0]), geom.pmf(x, geom_p, geom_pmin_low, geom_pmax_low))
        y_geom_high = np.minimum(np.ones(y_geom.shape[0]), geom.pmf(x, geom_p, geom_pmin_high, geom_pmax_high))

        fig = plt.figure(figsize=(12,8), dpi=100, facecolor='w', edgecolor='k')
        ax = fig.add_subplot(111)
        ax.xaxis.labelpad = 20
        ax.yaxis.labelpad = 20

        plt.plot(qlen, c2t,
                color='#bd0d45',
                alpha=0.1,
                label='C to T transitions')

        plt.plot(qlen, g2a,
                color='#236cf5',
                alpha=0.1,
                label='G to A transitions')

        plt.plot(x, y_unif, 
            linewidth=2.5, 
            color = 'DarkOliveGreen',
            alpha = 0.8,
            label = 'Uniform model')

        plt.fill_between(x, y_unif_low, y_unif_high,
            color='DarkOliveGreen',
            alpha=0.1,
            label = 'Uniform CI (2 sigma)')

        plt.plot(x, y_geom,
                linewidth=2.5, 
                color = '#D7880F',
                alpha = 0.8,
                label = 'Geometric model')

        plt.fill_between(x, y_geom_low, y_geom_high,
            color='#D7880F',
            alpha=0.1,
            label = 'Geometric CI (2 sigma)')

        plt.xlabel("Base from 5'", fontsize=20)
        plt.ylabel("Substitution frequency", fontsize=20)
        plt.suptitle(f"{contig}", fontsize=20, y=0.95)
        plt.title(f"coverage: {round(coverage,2)} | pvalue{rpval}",fontsize=12)
        plt.xticks(rotation=45, fontsize=8)
        ax.legend(fontsize=12)
        ax.set_xticks(qlen)
        plt.savefig(f"{plotdir}/{contig}.png", dpi=200)