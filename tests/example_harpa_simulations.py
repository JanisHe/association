"""
Example to simulate picks in an area and to pick these events with HARPA.
"""

import time

import numpy as np
from obspy import UTCDateTime
import pandas as pd
import seisbench.models as sbm  # noqa
import multiprocessing as mp

from simulation.arrivals import create_event, add_random_picks
from association.core.nll_functions import nll_wrapper
from association.core.interfaces import interface_harpa

from simulation.utils import moveout
import matplotlib.pyplot as plt


# Load stations
station_df = pd.read_csv("../metadata/stations.csv")

# Create picks dataframe
stations = pd.read_csv("../metadata/stations.csv")
lat_max = 48.90599442
lat_min = 48.89217758
lon_max = 7.938097954
lon_min = 7.914076328
max_depth = 3052.23012
min_depth = 1788.609028
starttime = UTCDateTime("20250401 15:30")

vp = 4500
vs = vp / np.sqrt(3)
true_events = ""

for i in range(10):  # Number of events
    origin_time = starttime + np.random.randint(0, 120)
    df, info = create_event(
        stations=stations,
        latitude_range=(lat_min, lat_max),
        longitude_range=(lon_min, lon_max),
        depth_range=(min_depth, max_depth),
        origin_time=origin_time,
        vp=vp,
        vs=vs,
    )
    if i == 0:
        df_events = df
    else:
        df_events = pd.concat(objs=[df_events, df])

    # Create string for true event information
    true_events += f"{origin_time} | {np.round(info[0], 3)}, {np.round(info[1], 3)}\n"

picks = add_random_picks(event_df=df_events, percentage=50)

# Start HARPA
verbose = 0

# Pick df for SeisBench catalogue data
pick_df = []
for idx in range(len(picks["trace_id"])):
    pick_df.append(
        {
            "id": picks.loc[idx, "trace_id"],
            "timestamp": UTCDateTime(picks.loc[idx, "peak_time"]).datetime,
            "prob": picks.loc[idx, "peak_value"],
            "type": picks.loc[idx, "phase"].lower(),
        }
    )
pick_df = pd.DataFrame(pick_df)

# Timing
stime = time.time()

# Start association
print("number of cpu:", mp.cpu_count())

config = {}
config["vel"] = {"P": vp / 1000, "S": vs / 1000}
config["ncpu"] = 8
config["P_phase"] = True
config["S_phase"] = True
config["z(km)"] = (0, 10)
config["min_peak_pre_event"] = 6
config["min_peak_pre_event_s"] = 0
config["min_peak_pre_event_p"] = 0
config["time_before"] = 5
config["DBSCAN"] = True
config["dbscan_min_samples"] = 1
config["epoch_before_decay"] = 2000  # Training more epochs can help find more events.
config["epoch_after_decay"] = 2000

stime = time.time()
catalog = interface_harpa(
    picks=pick_df, stations=station_df, config=config, verbose=verbose, iterations=1
)

print("Time needed for association:", time.time() - stime)

# Start localisation with NonLinLoc
catalog = nll_wrapper(
    catalog=catalog,
    station_df=station_df,
    nll_basepath="../NonLinLoc",
    nll_executable="/work/software/nlloc/7.00.16/src/bin",
    vel_model="../metadata/velocity_model.nll",
)

print("Time for association and NLL:", time.time() - stime)
print(catalog)

# Plot moveout curve
fig, ax = plt.subplots(nrows=1, ncols=1)
moveout(catalog=catalog, stations=stations, ax_lon=ax)
plt.show()
