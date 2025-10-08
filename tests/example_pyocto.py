import pandas as pd
import matplotlib.pyplot as plt

from association.core.interfaces import interface_pyocto
from association.core.nll_functions import nll_wrapper
from simulation.utils import moveout


# Load stations
station_df = pd.read_csv("../metadata/stations.csv")

# Read picks
pick_df = pd.read_csv("../metadata/picks.csv")

# Load velocity model
vel_model = pd.read_csv("../metadata/velocity_model.nll")

# Define config for PyOcto
config = {
    # "p_velocity": 4500,  # Not necessary when velocity model is used
    # "s_velocity": 2600,  # Not necessary when velocity model is used
    "zlim": (0, 30),
    "time_before": 10,
    "n_picks": 6,
    "n_p_picks": 3,
    "n_s_picks": 3,
    "n_p_and_s_picks": 3,
    "velocity_model_filename": "../metadata/velocity_model_pyocto",
}

# Start association
catalog = interface_pyocto(
    stations=station_df,
    picks=pick_df,
    config=config,
    velocity_model=vel_model,  # If p_- and s_velocity are defined in config set to None
    station_column_renaming={"trace_id": "id", "elevation_m": "elevation"},
    pick_column_renaming={"id": "station", "timestamp": "time", "type": "phase"},
)

# Start localisation with NonLinLoc
catalog = nll_wrapper(
    catalog=catalog,
    station_df=station_df,
    nll_basepath="../NonLinLoc",
    nll_executable="/work/software/nlloc/7.00.16/src/bin",
    vel_model="../metadata/velocity_model.nll",
)

print(catalog)

# Plot moveout curve
fig, ax = plt.subplots(nrows=1, ncols=1)
moveout(catalog=catalog, stations=station_df, ax_lon=ax)
plt.show()
