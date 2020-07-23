import pandas as pd
import numpy as np
import pyrolite.geochem
from pyrolite.util.pd import read_table
from pathlib import Path
import matplotlib.pyplot as plt
from pyrolite_meltsutil.tables.load import import_batch_config
from pyrolite_meltsutil.vis.style import COLORS

from pyrolite.util.plot.style import mappable_from_values
from pyrolite.util.plot.helpers import add_colorbar
import matplotlib.colors
from pyrolite.util.plot import save_figure, save_axes

outputfolder = Path("../data/experiments")

system, phases = (
    pd.read_hdf(outputfolder / "system.h5"),
    pd.read_hdf(outputfolder / "phases.h5"),
)
cfg = import_batch_config(outputfolder)

exprs = [
    hsh
    for (hsh, (name, c, e)) in cfg.items()
    if ((c["modes"] == ["isobaric"]) and (hsh in phases.experiment.unique()))
]

phases = phases.loc[phases.experiment.isin(exprs)]

ndexprs = phases.experiment.isin(
    [
        hsh
        for (hsh, (name, c, e)) in cfg.items()
        if (c["Title"].startswith("Nd") and (hsh in phases.experiment.unique()))
    ]
)

greenexprs = phases.experiment.isin(
    [
        hsh
        for (hsh, (name, c, e)) in cfg.items()
        if ((not c["Title"].startswith("Nd")) and (hsh in phases.experiment.unique()))
    ]
)

nexps = len(exprs)

exp_cpx = phases.phase.isin(["clinopyroxene"])
exp_opx = phases.phase.isin(["orthopyroxene"])

exp_spinel = phases.phase.isin(["spinel"])
exp_feldspar = phases.phase.isin(["feldspar"])
exp_olivine = phases.phase.isin(["olivine"])
exp_liquid = phases.phase.isin(["liquid"])
exp_cumulate = phases.phase.isin(["cumulate"])

#%% load real data
minchem = read_table("../data/all_traces.csv")
spinel_chem = read_table("../data/spinel_majors.csv")
spinel_chem.columns = [c.replace("_pct", "") for c in spinel_chem.columns]
minchem.columns = [c.replace("_pct", "") for c in minchem.columns]
minchem = minchem.rename(
    columns={"2Al/3O": "Al2O3", "2Cr/3O": "Cr2O3", "2V/5O": "V2O5"}
)
minchem.pyrochem.compositional = minchem.pyrochem.compositional.apply(
    pd.to_numeric, errors="coerce"
)
minchem = minchem.pyrochem.add_MgNo()
plag = minchem.Mineral == "plagioclase"
olivine = minchem.Mineral == "olivine"
clinopyroxene = minchem.Mineral == "pyroxene"
orthopyroxene = minchem.Mineral == "orthopyroxene"

#%%
title_kw = {"y": 0.95, "x": 0.3, "ha": "right", "va": "top"}
density_kw = dict(
    bins=100,
    contours=[0.9],
    label_contours=False,
    linestyles="--",
    linewidths=2,
    colors="r",
)
#%%
targets = ["Al2O3", "CaO", "MgO"]

fig, ax = plt.subplots(1, 2, figsize=(8, 4), subplot_kw={"projection": "ternary"})
ax[0].set_title("Clinopyroxene", **title_kw)
ax[1].set_title("Orthopyroxene", **title_kw)

for fltr, cmap in [(ndexprs, "plasma"), (greenexprs, "viridis")]:
    phases.loc[exp_cpx & fltr, targets].pyroplot.scatter(
        ax=ax[0], c=phases.loc[exp_cpx & fltr, "temperature"], s=3, cmap=cmap
    )
    phases.loc[exp_opx & fltr, targets].pyroplot.scatter(
        ax=ax[1], c=phases.loc[exp_opx & fltr, "temperature"], s=3, cmap=cmap
    )

minchem.loc[clinopyroxene, targets].pyroplot.density(ax=ax[0], **density_kw)
minchem.loc[orthopyroxene, targets].pyroplot.density(ax=ax[1], **density_kw)

plt.tight_layout()

save_figure(
    fig, name="Pyroxenes-Ternary", save_at="../img", save_fmts=["png", "pdf"], dpi=600
)
#%%
targets = ["Mg#", "Al2O3"]
fig, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(11, 4))
ax[0].set_title("Clinopyroxene", fontsize="x-large")
ax[1].set_title("Orthopyroxene", fontsize="x-large")
for fltr, cmap in [(ndexprs, "plasma"), (greenexprs, "viridis")]:
    phases.loc[exp_cpx & fltr, targets].pyroplot.scatter(
        ax=ax[0], c=phases.loc[exp_cpx & fltr, "temperature"], s=3, cmap=cmap
    )
    phases.loc[exp_opx & fltr, targets].pyroplot.scatter(
        ax=ax[1], c=phases.loc[exp_opx & fltr, "temperature"], s=3, cmap=cmap
    )

minchem.loc[clinopyroxene, targets].pyroplot.density(ax=ax[0], **density_kw)
minchem.loc[orthopyroxene, targets].pyroplot.density(ax=ax[1], **density_kw)
ax[0].set_ylim(0, 6)

plt.subplots_adjust(wspace=0.2)
cb = fig.colorbar(mappable_from_values(phases.loc[exp_opx, "temperature"]), ax=ax)
cb.set_label("Temperature", rotation=270, labelpad=20)
save_figure(fig, name="Pyroxenes", save_at="../img", save_fmts=["png", "pdf"], dpi=800)

#%%
targets = ["Al2O3", "CaO", "Na2O"]

fig, ax = plt.subplots(1)
ax.set_title("Plagoiclase", **title_kw)
ax = minchem.loc[plag, targets].pyroplot.density(ax=ax, **density_kw)
ax = phases.loc[exp_feldspar, targets].pyroplot.scatter(
    ax=ax, c=phases.loc[exp_feldspar, "temperature"],
)

save_figure(fig, name="Feldspar", save_at="../img", save_fmts=["png", "pdf"], dpi=800)
#%%
targets = ["Cr2O3", "Al2O3", "FeO"]

fig, ax = plt.subplots(1, 2, subplot_kw={"projection": "ternary"}, figsize=(8, 4))

for fltr, cmap in [(ndexprs, "plasma"), (greenexprs, "viridis")]:
    phases.loc[exp_spinel & fltr, targets].pyroplot.scatter(
        ax=ax[0], c=phases.loc[exp_spinel & fltr, "temperature"], s=3, cmap=cmap
    )

spinel_chem.loc[:, targets].pyroplot.density(ax=ax[0], **density_kw)

targets = ["Cr2O3", "Al2O3", "MgO"]
for fltr, cmap in [(ndexprs, "plasma"), (greenexprs, "viridis")]:
    phases.loc[exp_spinel & fltr, targets].pyroplot.scatter(
        ax=ax[1], c=phases.loc[exp_spinel & fltr, "temperature"], s=3, cmap=cmap
    )

spinel_chem.loc[:, targets].pyroplot.density(ax=ax[1], **density_kw)
plt.tight_layout()
save_figure(fig, name="Spinel", save_at="../img", save_fmts=["png", "pdf"], dpi=800)

#%%
targets = ["Mg#", "FeO"]
"""
ax = phases.loc[exp_liquid, targets].pyroplot.scatter(
    c=phases.loc[exp_liquid, "temperature"], cmap="plasma"
)
"""
ax = minchem.loc[olivine, targets].pyroplot.density(**density_kw)
phases.loc[exp_olivine, :].pyrochem.oxides.describe()
for fltr, cmap in [(ndexprs, "plasma"), (greenexprs, "viridis")]:
    phases.loc[exp_olivine & fltr, targets].pyroplot.scatter(
        ax=ax, c=phases.loc[exp_olivine & fltr, "temperature"], s=3, cmap=cmap
    )

ax.set_title("Olivine", **{**title_kw, "x": 0.9, "y": 0.9})

cb = add_colorbar(mappable_from_values(phases.loc[exp_olivine, "temperature"]), ax=ax,)
cb.set_label("Temperature", rotation=270, labelpad=20)

save_figure(
    ax.figure, name="Olivine", save_at="../img", save_fmts=["png", "pdf"], dpi=600
)
ax.set(xlim=(0.65, 0.9), ylim=(14, 27))
save_figure(
    ax.figure, name="Olivine-Zoom", save_at="../img", save_fmts=["png", "pdf"], dpi=600
)
#%%
targets = ["Al2O3", "MgO", "FeO"]
fig, ax = plt.subplots(2, 1, figsize=(4, 8), subplot_kw={"projection": "ternary"})

ax[0].set_title("Liqiud", **title_kw)
ax[1].set_title("Cumulate", **title_kw)

for fltr, cmap in [(ndexprs, "plasma"), (greenexprs, "viridis")]:
    phases.loc[exp_liquid & fltr, targets].pyroplot.scatter(
        ax=ax[0], c=phases.loc[exp_liquid & fltr, "temperature"], s=3, cmap=cmap
    )
    phases.loc[exp_cumulate & fltr, targets].pyroplot.scatter(
        ax=ax[1], c=phases.loc[exp_cumulate & fltr, "temperature"], s=3, cmap=cmap
    )
plt.tight_layout()

save_figure(
    fig, name="Liquid-Cumulate", save_at="../img", save_fmts=["png", "pdf"], dpi=600
)
