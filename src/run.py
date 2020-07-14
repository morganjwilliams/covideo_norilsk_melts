from pathlib import Path
from pyrolite.util.meta import stream_log
from pyrolite.util.pd import read_table

logger = stream_log("pyrolite-meltsutil", level="INFO")
outputfolder = Path("../data/experiments")

df = read_table("../data/starting_compositions.csv").rename(
    columns={"FeOt": "FeO", "LOI": "H2O"}
)

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
