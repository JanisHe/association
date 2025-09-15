import pandas as pd

from core.interfaces import interface_pyocto
from core.nll_functions import nll_wrapper


# Load stations
station_df = pd.read_csv("../metadata/stations.csv")

# Read picks
pick_df = pd.read_csv("../metadata/picks.csv")

# Define config for PyOcto
config = {
    "p_velocity": 4500,
    "s_velocity": 2600,
    "zlim": (0, 30),
    "time_before": 10,
    "n_picks": 6,
    "n_p_picks": 3,
    "n_s_picks": 3,
    "n_p_and_s_picks": 3,
}

# Start association
catalog = interface_pyocto(
    stations=station_df,
    picks=pick_df,
    velocity_model=None,
    **config
)

# Start localisation with NonLinLoc
catalog = nll_wrapper(catalog=catalog,
                      station_df=station_df,
                      nll_basepath="../NonLinLoc",
                      nll_executable="/work/software/nlloc/7.00.16/src/bin",
                      vel_model="../metadata/velocity_model.nll")

print(catalog)
