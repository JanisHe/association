import pandas as pd

from core.interfaces import interface_harpa
from core.nll_functions import nll_wrapper


# Load stations
station_df = pd.read_csv("../metadata/stations.csv")

# Read picks
pick_df = pd.read_csv("../metadata/picks.csv")

config = {
    "vel": {"P": 6.0, "S": 6.0 / 1.75},
    "ncpu": 6,
    "P_phase": True,
    "S_phase": True,
    "z(km)": (0, 10),
    "min_peak_pre_event": 6,
    "min_peak_pre_event_p": 0,
    "min_peak_pre_event_s": 0,
    "time_before": 5,
    "DBSCAN": True,
    "dbscan_min_samples": 1,
    "epoch_before_decay": 2000,  # Training more epochs can help find more events
    "epoch_after_decay": 2000
}

catalog = interface_harpa(
    picks=pick_df,
    stations=station_df,
    config=config,
    verbose=0,
    iterations=1
)

# Start localisation with NonLinLoc
catalog = nll_wrapper(catalog=catalog,
                      station_df=station_df,
                      nll_basepath="../NonLinLoc",
                      nll_executable="/work/software/nlloc/7.00.16/src/bin",
                      vel_model="../metadata/velocity_model.nll")

print(catalog)
