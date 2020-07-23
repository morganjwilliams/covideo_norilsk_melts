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
        ((c["modes"] == ["isobaric"]) and c["Title"] == "Mr2")
        and (hsh in phases.experiment.unique())
    )
]

exprs = sorted(exprs, key=lambda e: cfg[e][1]["H2O"])  # sort by H2O content

fig, ax = _phasevolumes(
    phases.loc[phases.experiment.isin(exprs)],
    config=cfg,
    n_across=3,
    aspect=1,
    unit_size=5,
    exprs=exprs,
)
for a, expr in zip(ax, exprs):
    name, config, env = cfg[expr]
    a.set_title(
        "$\mathrm{H_2O="
        + "{:.1f}".format(np.round(config.get("H2O", 0), 1))
        + "}$ Wt%",
        y=1.05,
        fontsize="x-large",
    )

save_figure(
    fig, name="Mr2_H2O_variation", save_at="../img/", save_fmts=["png", "pdf"], dpi=800
)
