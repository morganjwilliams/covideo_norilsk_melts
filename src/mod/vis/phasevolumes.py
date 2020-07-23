import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from pyrolite.util.plot.legend import proxy_line
from pyrolite_meltsutil.vis.style import phase_color, phaseID_linestyle
from pyrolite_meltsutil.vis.templates import plot_phasevolumes
from pyrolite_meltsutil.util.tables import phasename
from ..sequence import get_appearance_sequence


def _phasevolumes(
    phases,
    config=None,
    n_across=4,
    tlim=(800, 1250),
    unit_size=5,
    aspect=0.8,
    exprs=None,
    legend_on=None,
):
    if exprs is None:
        exprs = phases["experiment"].unique()
    ndown = len(exprs) // n_across
    if len(exprs) % n_across:
        ndown += 1
    figsize = (n_across * unit_size, aspect * ndown * unit_size)
    fig, ax = plt.subplots(ndown, n_across, sharex=True, sharey=True, figsize=figsize,)
    ax = ax.flat
    ax[0].set_yscale("log")
    ax[0].set_ylim((0.1, 100))
    ax[0].set_xlim(tlim)
    sequences = []
    for ix, expnumber in enumerate(exprs):
        name = expnumber
        exp_title, exp_config, exp_env = config[name]
        expdf = phases.loc[phases.experiment == expnumber, :]
        sequence = get_appearance_sequence(expdf)
        sequences.append(
            {
                phs: sequence.index((v, _phases))
                for (v, _phases) in sequence
                for phs in _phases
            }
        )
        plot_phasevolumes(expdf, marker=None, ax=ax[ix], lw=4, legend=False)
        ax[ix].set_title(exp_config["Suite"] + ": " + exp_config["Title"])

        modes = exp_config["modes"]
    all_phaseIDs = sorted(
        [
            (p, np.mean([s.get(p, len(s.keys()) + 1) for s in sequences]))
            for p in phases.phaseID.unique()
            if not pd.isnull(p)
        ],
        key=lambda x: x[1],
    )
    # convert IDs to names where its the first or only
    name_counter = Counter([phasename(pID) for (pID, pos) in all_phaseIDs])
    all_phaseIDs = [
        (phasename(pID), pos)
        if (name_counter[phasename(pID)] == 1) or (pID.endswith("_0"))
        else (pID, pos)
        for (pID, pos) in all_phaseIDs
    ]
    proxies = {
        k: proxy_line(color=phase_color(k), ls=phaseID_linestyle(k), lw=4)
        for k, pos in all_phaseIDs
    }
    if legend_on is None:
        ax[n_across - 1].legend(proxies.values(), proxies.keys(), fontsize="large")
    else:
        ax[legend_on].legend(proxies.values(), proxies.keys(), fontsize="large")
    plt.subplots_adjust(hspace=0.3)
    return fig, ax


"""
ax[ix].annotate(modes, xy=(0.1, 0.9), xycoords="axes fraction", ha="left")

h2o = exp_config["H2O"]
ax[ix].annotate(
    "H2O={}".format(h2o), xy=(0.1, 0.8), xycoords="axes fraction", ha="left"
)
"""
