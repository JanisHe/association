"""
Compare different associators
"""

import time
import pandas as pd
import matplotlib.pyplot as plt

from association.core.interfaces import (
    interface_pyocto,
    interface_gamma,
    interface_harpa,
)
from simulation.utils import moveout


# Load stations
station_df = pd.read_csv("../metadata/stations.csv")

# Read picks
pick_df = pd.read_csv("../metadata/picks.csv")

iterations = 3  # Number of second pass iterations

# Config PyOcto
config_pyocto = {
    "p_velocity": 4500,  # Not necessary when velocity model is used
    "s_velocity": 2600,  # Not necessary when velocity model is used
    "zlim": (0, 30),
    "time_before": 10,
    "n_picks": 6,
    "n_p_picks": 3,
    "n_s_picks": 3,
    "n_p_and_s_picks": 3,
    "second_pass_overwrites": {
        "time_before": 10,
        "n_picks": 6,
        "n_p_picks": 3,
        "n_s_picks": 3,
        "n_p_and_s_picks": 3,
        "iterations": iterations,
    },
}

# Config GaMMA
config_gamma = {
    "vel": {"p": 4.5, "s": 2.6},
    "z(km)": (0, 30),
    "ncpu": 8,
    "use_dbscan": True,
    "use_amplitude": False,
    "min_picks_per_eq": 6,
    "min_p_picks_per_eq": 3,
    "min_s_picks_per_eq": 3,
    "max_sigma11": 2.0,
}


# Config HARPA
config_harpa = {
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
    "epoch_after_decay": 2000,
}

# Build catalog for each associator
stime_harpa = time.time()
catalog_harpa = interface_harpa(
    picks=pick_df,
    stations=station_df,
    config=config_harpa,
    verbose=0,
    second_pass_iterations=iterations,
)
time_harpa = time.time() - stime_harpa

stime_pyocto = time.time()
catalog_pyocto = interface_pyocto(
    stations=station_df,
    picks=pick_df,
    velocity_model=None,
    config=config_pyocto,
    station_column_renaming={"trace_id": "id", "elevation_m": "elevation"},
    pick_column_renaming={"id": "station", "timestamp": "time", "type": "phase"},
)
time_pyocto = time.time() - stime_pyocto

stime_gamma = time.time()
catalog_gamma = interface_gamma(
    stations=station_df,
    picks=pick_df,
    config=config_gamma,
    second_pass_iterations=iterations,
)
time_gamma = time.time() - stime_gamma

# Print number of events of each derived catalog
print("Events PyOcto:", len(catalog_pyocto), "| Time PyOcto:", time_pyocto)
print("Events GaMMA:", len(catalog_gamma), "| Time GaMMA:", time_gamma)
print("Events HARPA:", len(catalog_harpa), "| Time HARPA:", time_harpa)

# Plot moveout curve for all catalogs
fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=True)
moveout(catalog=catalog_pyocto, stations=station_df, ax_lon=ax[0], title="PyOcto")
moveout(catalog=catalog_gamma, stations=station_df, ax_lon=ax[1], title="GaMMA")
moveout(catalog=catalog_harpa, stations=station_df, ax_lon=ax[2], title="HARPA")
plt.show()
