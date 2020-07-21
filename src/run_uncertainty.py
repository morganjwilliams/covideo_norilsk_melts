from pathlib import Path
import numpy as np
from pyrolite.util.meta import stream_log
from pyrolite.util.pd import read_table
from pyrolite.comp.codata import ILR, inverse_ILR
import pyrolite.geochem
from pyrolite.util.text import slugify
from pyrolite.util.pd import accumulate
from pyrolite.util.plot import save_figure

np.random.seed(32)


def blur_compositions(df, noise=0.05, scale=100):
    """
    Function to add 'compositional noise' to a set of compositions. In reality, it's
    its best to use measured uncertainties to generate these simulated compositions.
    """
    # transform into compositional space, add noise, return to simplex
    xvals = ILR(df.values)
    xvals += np.random.randn(*xvals.shape) * noise
    return inverse_ILR(xvals) * scale


logger = stream_log("pyrolite-meltsutil", level="INFO")
outputfolder = Path("../data/experiments_uncertainty")

# Naldret 2004 compositions
df = read_table("../data/starting_compositions.csv").rename(
    columns={"FeOt": "FeO", "LOI": "H2O"}
)

reps = 10  # increase this to perform more experiments
df = accumulate([df.iloc[[-1], :]] * reps)
df = df.reset_index()  # .drop(columns="index")
df.pyrochem.compositional = blur_compositions(df.pyrochem.compositional)
ax = df[["Al2O3", "CaO", "SiO2"]].pyroplot.scatter(s=5, c="k", alpha=0.5)
save_figure(
    ax.figure,
    name="Mr2_Uncertainty_Distribution",
    save_at="../img/",
    save_fmts=["png", "pdf"],
)
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
