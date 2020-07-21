import pandas as pd
import numpy as np
import pyrolite.geochem
from pyrolite.util.pd import read_table
from pathlib import Path
import matplotlib.pyplot as plt
from pyrolite_meltsutil.tables.load import import_batch_config
from pyrolite_meltsutil.vis.style import COLORS

import matplotlib.colors
from pyrolite.util.plot import save_figure

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
phases.loc[exp_cpx, targets].pyroplot.scatter(
    ax=ax[0], c=phases.loc[exp_cpx, "temperature"], cmap="Greens", s=3
)
phases.loc[exp_opx, targets].pyroplot.scatter(
    ax=ax[1], c=phases.loc[exp_opx, "temperature"], cmap="Purples", s=2
)
minchem.loc[clinopyroxene, targets].pyroplot.density(ax=ax[0], **density_kw)
minchem.loc[orthopyroxene, targets].pyroplot.density(ax=ax[1], **density_kw)

plt.tight_layout()

save_figure(fig, name="Pyroxenes", save_at="../img", save_fmts=["png", "pdf"])
#%%
targets = ["Al2O3", "CaO", "Na2O"]

fig, ax = plt.subplots(1)
ax.set_title("Plagoiclase", **title_kw)
ax = minchem.loc[plag, targets].pyroplot.density(ax=ax, **density_kw)
ax = phases.loc[exp_feldspar, targets].pyroplot.scatter(
    ax=ax, c=phases.loc[exp_feldspar, "temperature"],
)

save_figure(fig, name="Feldspar", save_at="../img", save_fmts=["png", "pdf"])
#%%
targets = [
    "SiO2",
    "MgO",
    "FeO",
]

fig, ax = plt.subplots(1)

ax = minchem.loc[olivine, targets].pyroplot.density(ax=ax, **density_kw)
ax = phases.loc[exp_olivine, targets].pyroplot.scatter(
    ax=ax, c=phases.loc[exp_olivine, "temperature"]
)
ax.set_title("Olivine", **title_kw)
save_figure(fig, name="Olivine", save_at="../img", save_fmts=["png", "pdf"])
#%%
targets = ["Al2O3", "TiO2", "FeO"]

fig, ax = plt.subplots(1)

ax = phases.loc[exp_spinel, targets].pyroplot.scatter(
    ax=ax, c=phases.loc[exp_spinel, "temperature"]
)
spinel_chem.loc[:, targets].pyroplot.density(ax=ax, **density_kw)

ax.set_title("Spinel", **title_kw)
save_figure(fig, name="Spinel", save_at="../img", save_fmts=["png", "pdf"])
phases.loc[exp_liquid, targets].pyroplot.scatter(
    ax=ax, c=phases.loc[exp_liquid, "temperature"], cmap="plasma"
)
ax.set_title("Spinel+Liquid", **title_kw)
save_figure(fig, name="Spinel+Liquid", save_at="../img", save_fmts=["png", "pdf"])
#%%
from pyrolite.util.plot.style import mappable_from_values
from pyrolite.util.plot.helpers import add_colorbar

targets = ["MgO", "FeO"]
ax = phases.loc[exp_liquid, targets].pyroplot.scatter(
    c=phases.loc[exp_liquid, "temperature"], cmap="plasma"
)
phases.loc[exp_olivine, targets].pyroplot.scatter(
    c=phases.loc[exp_olivine, "temperature"], ax=ax
)
minchem.loc[olivine, targets].pyroplot.density(ax=ax, **density_kw)
ax.set_title("Olivine+Liquid", **{**title_kw, "x": 0.9, "y": 0.9})

add_colorbar(
    mappable_from_values(phases.loc[exp_olivine, "temperature"]), ax=ax,
)
"""
add_colorbar(
    mappable_from_values(phases.loc[exp_liquid, "temperature"], cmap="plasma"), ax=ax,
)
"""
save_figure(
    ax.figure, name="Liquid-Olivine", save_at="../img", save_fmts=["png", "pdf"]
)
#%%
targets = ["Al2O3", "MgO", "FeO"]
fig, ax = plt.subplots(1, 2, figsize=(8, 4), subplot_kw={"projection": "ternary"})
ax[0].set_title("Cumulate", **title_kw)
ax[1].set_title("Liqiud", **title_kw)
phases.loc[exp_cumulate, targets].pyroplot.scatter(
    ax=ax[0], c=phases.loc[exp_cumulate, "temperature"],
)
phases.loc[exp_liquid, targets].pyroplot.scatter(
    ax=ax[1], c=phases.loc[exp_liquid, "temperature"], cmap="plasma"
)
plt.tight_layout()

save_figure(fig, name="Liquid-Cumulate", save_at="../img", save_fmts=["png", "pdf"])
