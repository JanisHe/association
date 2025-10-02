import pandas as pd
import matplotlib.pyplot as plt

from association.core.interfaces import interface_gamma
from association.core.nll_functions import nll_wrapper
from simulation.utils import moveout


# Load stations
station_df = pd.read_csv("../metadata/stations.csv")

# Read picks
pick_df = pd.read_csv("../metadata/picks.csv")

# Read velocity model
velocity_model = pd.read_csv("../metadata/velocity_model.nll")

# Define config for GaMMA
# Note more parameters are available (https://ai4eps.github.io/GaMMA/)
config = {
    "vel": {
        "p": 4.5,  # If velocity model is used then the velocity model will be used
        "s": 2.6,
    },
    "z(km)": (0, 30),
    "ncpu": 4,
    "use_dbscan": True,
    "use_amplitude": False,
    "min_picks_per_eq": 6,
    "min_p_picks_per_eq": 3,
    "min_s_picks_per_eq": 3,
    "max_sigma11": 2.0,
}

# Start association
catalog = interface_gamma(
    stations=station_df,
    picks=pick_df,
    velocity_model=velocity_model,
    config=config,
    verbose=True,
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
