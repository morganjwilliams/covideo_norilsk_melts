import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

import pyrolite.geochem
from pyrolite.util.pd import read_table

from pyrolite.util.plot import save_figure
from pyrolite_meltsutil.tables.load import import_batch_config
from pyrolite_meltsutil.vis.style import COLORS
from mod.vis.phasevolumes import _phasevolumes

outputfolder = Path("../data/experiments")

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
        (c["modes"] == ["isobaric"] and c["H2O"] == 0)
        and (hsh in phases.experiment.unique())
    )
]

fig, ax = _phasevolumes(
    phases.loc[phases.experiment.isin(exprs)],
    config=cfg,
    n_across=4,
    aspect=1,
    unit_size=5,
)
save_figure(fig, name="Batch_Dry", save_at="../img/", save_fmts=["png", "pdf"])
#%%
exprs = [
    hsh
    for (hsh, (name, c, e)) in cfg.items()
    if (
        (c["modes"] == ["isobaric"] and c["H2O"] == 1)
        and (hsh in phases.experiment.unique())
    )
]
fig, ax = _phasevolumes(
    phases.loc[phases.experiment.isin(exprs)],
    config=cfg,
    n_across=4,
    aspect=1,
    unit_size=5,
)
save_figure(fig, name="Batch_1Wt%H2O", save_at="../img/", save_fmts=["png", "pdf"])

#%%

exprs = [
    hsh
    for (hsh, (name, c, e)) in cfg.items()
    if (
        (c["modes"] == ["isobaric", "fractionate solids"] and c["H2O"] == 1)
        and (hsh in phases.experiment.unique())
    )
]
fig, ax = _phasevolumes(
    phases.loc[phases.experiment.isin(exprs)],
    config=cfg,
    n_across=4,
    aspect=1,
    unit_size=5,
)
save_figure(fig, name="Frac_1wt%H2O", save_at="../img/", save_fmts=["png", "pdf"])

#%%
#%%
from mod.sequence import (
    get_appearance_sequence,
    get_assemblage_sequence,
    sequence_distance,
)

sequence_distance(
    phases.loc[phases.experiment == exprs[0]],
    ["liquid", "clinopyroxene", "feldspar", "olivine"],
    ignore_trailing=True,
)

{
    e: get_appearance_sequence(phases.loc[phases.experiment == e])
    for e in phases.experiment.unique()
}

# could use graph edit distance for measuring distance between sequences where
# multiple things can come in at once

better_models = sorted(
    [
        (
            e,
            sequence_distance(
                phases.loc[phases.experiment == e],
                ["liquid", "olivine", "clinopyroxene", "orthopyroxene", "feldspar",],
                ignore_trailing=True,
            )[-1],
        )
        for e in phases.experiment.unique()
    ],
    key=lambda x: x[1],
)

fig, ax = plt.subplots(1)

ax.hist2d(*np.array([(dist, cfg[e][1]["H2O"]) for (e, dist) in better_models]).T,)
