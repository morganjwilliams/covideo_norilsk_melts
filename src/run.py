from pathlib import Path
import numpy as np
from pyrolite.util.meta import stream_log
from pyrolite.util.pd import read_table
import matplotlib.pyplot as plt
from pyrolite.util.plot import save_figure
from pyrolite.util.plot.legend import proxy_line
from pyrolite.util.units import scale

logger = stream_log("pyrolite-meltsutil", level="INFO")
outputfolder = Path("../data/experiments")

# Naldrett 2004 compositions
df = read_table("../data/starting_compositions_all.csv").rename(
    columns={"FeOt": "FeO", "LOI": "H2O"}
)
df.pyrochem.elements *= scale("ppm", "wt%")
df = df.pyrochem.convert_chemistry(
    [i for i in df.pyrochem.oxides if "Fe" not in i]
    + [{"FeO": 0.9, "Fe2O3": 0.1}]
    + ["Cr2O3"],
    renorm=True,
)

targets = ["Al2O3", "CaO", "SiO2"]

kw = dict(s=30, edgecolors="k")
fig, ax = plt.subplots(1, subplot_kw={"projection": "ternary"})
labels, proxies = [], []
styles = {
    "Nadezhdinsky": {"c": "yellow", "marker": "D"},
    "Morongovsky": {"c": "lime", "marker": "D"},
    "Mokulaevsky": {"c": "lime", "marker": "s"},
    "Kharaelakhsky": {"c": "lime", "marker": "o"},
}

for s in df.Suite.unique():
    ax = df.loc[df.Suite == s, targets].pyroplot.scatter(ax=ax, **{**kw, **styles[s]})
    labels.append(s)
    proxies.append(proxy_line(markersize=np.sqrt(kw["s"]), ls="", mec="k", **styles[s]))
ax.legend(proxies, labels, fontsize="x-large", markerscale=2, bbox_to_anchor=(0.75, 1))
save_figure(fig, name="StartingPoints", save_at="../img/", save_fmts=["png", "pdf"], dpi=800)

#%%
from pyrolite_meltsutil.env import MELTS_Env
from pyrolite_meltsutil.automation import MeltsBatch

env = MELTS_Env()
env.VERSION = "MELTS"
env.MODE = "isobaric"
env.MINP = 100
env.MAXP = 10000
env.MINT = 500
env.MAXT = 1500
env.DELTAT = -10
env.DELTAP = 0
env.spec["MAXP"]
batch = MeltsBatch(
    df,
    default_config={
        "Initial Pressure": 500,
        "Initial Temperature": 1250,
        "Final Temperature": 800,
        "modes": ["isobaric", "fractionate solids"],
    },
    config_grid={
        "Log fO2 Path": ["NNO"],
        "modes": [None, ["isobaric"]],
        "modifychem": [None, {"H2O": 0}, {"H2O": 1}],
    },
    env=env,
    fromdir=outputfolder,
    logger=logger,
)

batch.run(
    overwrite=False
)  # overwrite=False if you don't want to update existing exp folders
