import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

import pyrolite.geochem
from pyrolite.util.pd import read_table

from pyrolite.util.plot import save_figure
from pyrolite_meltsutil.tables.load import import_batch_config
from pyrolite_meltsutil.vis.style import (
    COLORS,
    phaseID_linestyle,
    phaseID_marker,
    phase_color,
)
from pyrolite.util.plot.legend import proxy_line

from mod.vis.phasevolumes import _phasevolumes
from mod.sequence import (
    get_appearance_sequence,
    get_assemblage_sequence,
    sequence_distance,
)
from pyrolite_meltsutil.util.tables import phasename
from collections import Counter

outputfolder = Path("../data/experiments_uncertainty")

system, phases = (
    pd.read_hdf(outputfolder / "system.h5"),
    pd.read_hdf(outputfolder / "phases.h5"),
)
cfg = import_batch_config(outputfolder)
#%%
exprs = [
    hsh
    for (hsh, (name, c, e)) in cfg.items()
    if (
        (c["modes"] == ["isobaric"] and c["H2O"] == 1)
        and (hsh in phases.experiment.unique())
    )
]

phases.phase.unique()
phaselist = [
    "liquid",
    "olivine",
    "clinopyroxene",
    "feldspar",
    "orthopyroxene",
    "spinel",
]

fig, ax = plt.subplots(
    len(phaselist) // 2, 2, sharex=True, sharey=True, figsize=(8, 8)
)
xvar, yvar = "temperature", "volume%"
[a.set_xlabel(xvar) for a in ax[-1, :]]
[a.set_ylabel(yvar) for a in ax[:, 0]]
name_counter = Counter(
    [
        phasename(pID)
        for pID in phases.phaseID.unique()[~pd.isnull(phases.phaseID.unique())]
    ]
)
for p, pax in zip(phaselist, ax.flat):
    pdf = phases.loc[(phases.phase == p) & (phases.experiment.isin(exprs)), :]
    proxies = {}
    for phaseID in pdf.phaseID.unique():
        style = dict(ls=phaseID_linestyle(phaseID), color=phase_color(phaseID))
        for expr in pdf.experiment.unique():
            e_p_df = pdf.loc[((pdf.phaseID == phaseID) & (pdf.experiment == expr)), :]
            pax.plot(e_p_df[xvar], e_p_df[yvar], **style)
            name = (
                phasename(phaseID)
                if ((name_counter[phasename(phaseID)] == 1) or (phaseID.endswith("_0")))
                else phaseID
            )
            proxies[name] = proxy_line(**style)

    pax.legend(
        proxies.values(),
        proxies.keys(),
        bbox_to_anchor=None,
        loc="best",
        fontsize="large",
    )

plt.tight_layout()

save_figure(
    fig, name="Batch_Uncertainty_1Wt%H2O", save_at="../img/", save_fmts=["png", "pdf"]
)
#%%

better_models = sorted(
    [
        (
            e,
            sequence_distance(
                phases.loc[phases.experiment == e],
                ["liquid", "olivine", "feldspar", "clinopyroxene", "orthopyroxene",],
                ignore_trailing=True,
            )[-1],
        )
        for e in phases.experiment.unique()
    ],
    key=lambda x: x[1],
)

fig, ax = plt.subplots(1)

ax.hist2d(*np.array([(dist, cfg[e][1]["H2O"]) for (e, dist) in better_models]).T,)
