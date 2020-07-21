from pathlib import Path
from pyrolite_meltsutil.tables.load import (
    aggregate_tables,
    import_batch_config,
    import_tables,
)
from pyrolite.util.meta import stream_log
stream_log('pyrolite_meltsutil')

outputfolder = Path("../data/experiments")
system, phases = aggregate_tables(outputfolder)

system.to_hdf(outputfolder / "system.h5", key="df")
phases.to_hdf(outputfolder / "phases.h5", key="df")

#%%
outputfolder = Path("../data/experiments_uncertainty")
system, phases = aggregate_tables(outputfolder)

system.to_hdf(outputfolder / "system.h5", key="df")
phases.to_hdf(outputfolder / "phases.h5", key="df")
