# Seismic Phase Association with different Associators

This packages contains scripts and interfaces to run either [PyOcto](https://github.com/yetinam/pyocto),
[HARPA](https://github.com/DaDaCheng/phase_association/tree/main) or
[GaMMA](https://github.com/AI4EPS/GaMMA) for seismic phase association on seismic picks. Seismic data
can be picked, for example, by using [PhaseNet](https://github.com/AI4EPS/PhaseNet)
or available Deep-Learning pickers from [SeisBench](https://github.com/seisbench/seisbench).

To learn how to pick seismic data, either follow the examples given in
[SeisBench](https://github.com/seisbench/seisbench) or use
[my package](https://github.com/JanisHe/seisbench_picking)
to pick seismic phase onsets with SeisBench.

## Required packages
* standard packages: `numpy`, `obspy`, `pandas`, `seisbench`, `matplotlib`
* Packages for phase association:
  - `PyOcto`: `pip install pyocto`,
     PyOcto requires Pyrocko: `https://github.com/pyrocko/pyrocko`
  - `HARPA`: `pip install -q git+https://github.com/DaDaCheng/phase_association.git`, <br>
     HARPA requires POT: `pip install POT` <br>
     HARPA requires scikit-learn: `pip install scikit-learn`
  - `GaMMA`: `pip install git+https://github.com/wayneweiqiang/GaMMA.git`

## Preparing files
### Stations
A `.csv` file that contains station information (note the file header!):

| id         | latitude | longitude | elevation |
|------------|----------|-----------|-----------|
| FO.BETS.00 | 48.89357 | 7.92429   | 146       |
| RG.RITT.00 | 48.89436 | 7.96103   | 138       |
| RG.KUHL.00 | 48.91473 | 7.92996   | 176       |
| FO.OPS.00  | 48.92126 | 7.88278   | 198       |

- `id`: ID of each station (`network.station.location`)\
   If `location` is not known use `network.station`
- `latitude`: Latitude of station in degree
- `longitude`: Longitude of station in degree
- `elevation`: Elevation of station in m above sea level

Since the column for the table above are only valid for HARPA, the column
names will be changed for PyOcto. This is done by defining a dictionary
as an argument of the function `interface_pytocto` (see example below).

### PhaseNet picks
A `.csv` file that contains all picks from all available stations:

| id         | timestamp                  | prob    | type |
|------------|----------------------------|---------|------|
| FO.BETS.00 | 2025-04-01 15:31:41.065029 | 0.48751 | P    |
| RG.RITT.00 | 2025-04-01 15:31:41.12786  | 0.76287 | P    |
| RG.RITT.00 | 2025-04-01 15:31:41.496281 | 0.62723 | S    |
| RG.KUHL.00 | 2025-04-01 15:31:42.015891 | 0.89125 | P    |
| RG.KUHL.00 | 2025-04-01 15:31:42.158434 | 0.68523 | S    |

- `id`: ID of station (`network.station.location`). Must fit with IDs of station file. \
   If `location` is not known use `network.station`
- `timestamp`: Time of pick from PhaseNet
- `prob`: Pick probability of PhaseNet
- `type`: Phase type (Either P or S)

For PyOcto, the column names will be changed. This is done by defining
a dictionary as an argument of the function `interface_pytocto` (see
example below).

### Velocity Model (Optional for PyOcto, required for NonLinLoc)
A `.csv` file that contains a velocity model of the area. If the velocity model is unknown,
define an own velocity model. A detailed velocity model is not required for this purpose.
Otherwise, the velocity model only has constant P- and S-velocities, which is also the case
for HARPA.

| depth  | vp    | vs    |
|--------|-------|-------|
| -0.55  | 1.318 | 0.621 |
| 0.074  | 2.279 | 1.073 |
| 0.418  | 2.551 | 1.202 |
| 0.65   | 2.954 | 1.52  |

- `depth`: Depth in km (positive downwards)
- `vp`: P-velocity in km/s
- `vs`: S-velocity in km/s

## Running Phase Association
Once all files are created, read the different `.csv` files by using `pandas'` function `read_csv`.
Then, phase association is started by running the following codes for the different associators
(config dictionaries can be varied by user; use available keywords that can be found in
the descriptions of [Pyocto](https://github.com/yetinam/pyocto), [Harpa](https://github.com/DaDaCheng/phase_association/tree/main) or [GaMMA](https://github.com/AI4EPS/GaMMA)):

### PyOcto
In comparison to HARPA and GaMMA, PyOcto needs other keywords. Therefore, the columns of the given
`.csv` files for the stations and picks are modified by using the arguments `station_colum_renaming`
and `pick_column_renaming`. The key of the dictionary represents the original column name and the
value is the required column name for PyOcto.

```python
from association.core.interfaces import interface_pyocto

# Association with Pyocto
config = {
    "p_velocity": vp,
    "s_velocity": vs,
    "zlim": (0, 30),
    "time_before": 10,
    "n_picks": 6,
    "n_p_picks": 3,
    "n_s_picks": 3,
    "n_p_and_s_picks": 3,
}

catalog = interface_pyocto(
    picks=picks,
    stations=stations,
    config=config,
    velocity_model=None,  # i.e. no velocity model is given
    station_column_renaming={"trace_id": "id", "elevation_m": "elevation"},
    pick_column_renaming={"id": "station", "timestamp": "time", "type": "phase"}
)
```

### HARPA
```python
from association.core.interfaces import interface_harpa


# Association with HARPA
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
    picks=picks,
    stations=stations,
    config=config
)
```

### GaMMA
```python
from association.core.interfaces import interface_gamma

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

catalog = interface_gamma(
    stations=station_df,
    picks=pick_df,
    velocity_model=velocity_model,
    config=config,
    verbose=True,
)
```

All interface functions return an `obspy.Catalog` that contains all events, including
the picks.

### Associating Remaining Picks
After the first association some picks remain. In order to associate events from these
remaining picks, a `second_pass_iteration` is possible for each associator. This is done
in PyOcto by modifying the config dictionary, for example:
```python
config_pyocto = {
    "p_velocity": 4500,  # Not necessary when velocity model is used
    "s_velocity": 2600,  # Not necessary when velocity model is used
    "zlim": (0, 30),
    "time_before": 10,
    "n_picks": 6,
    "n_p_picks": 3,
    "n_s_picks": 3,
    "n_p_and_s_picks": 3,
    "second_pass_overwrites": {  # Do a second or third association with remaining picks
        "time_before": 10,
        "n_picks": 6,
        "n_p_picks": 3,
        "n_s_picks": 3,
        "n_p_and_s_picks": 3,
        "iterations": 1,  # Number of iterations for association of remaining picks
    },
}
```
In HARPA and GaMMA both interfaces have the argument `second_pass_iterations` which is an
`int` and controls the number of additional iterations to associate the remaining picks.

## Relocalisation with NonLinLoc
After seismic phase association, the events can be relocalised using [NonLinLoc](http://alomax.free.fr/nlloc/).
Therefore, a velocity model (see above) is required.
In the following, two different approaches are described to use NLL
- semi-automatic approach: does relocalisation using a single function
- manual approach: creating `.obs` file(s) and localise with own methods

### Semi-Automatic Approach
```python
from association.core.nll_functions import nll_wrapper

catalog = nll_wrapper(catalog=catalog,
                      station_df=stations,
                      nll_basepath="../NonLinLoc",
                      nll_executable="/work/software/nlloc/7.00.16/src/bin",
                      vel_model="../metadata/rittershoffen.nll")
```

- `catalog`: Seismicity catalog from association
- `station_df`: pandas Dataframe with station information (same as above)
- `nll_basepath`: Pathname were all NLL files are saved
- `nll_executable`: Pathname of executable NLL files
- `vel_model`: Pathname (to .csv-file) or pandas Dataframe of velocity model

If `nll_basepath` already exists, the `time` and `model` files will not be
created by NonLinLoc. This is only done during the first run, when `nll_basepath`
is created by `nll_wrapper`.

### Manual Approach
```python
from association.core.nll_functions import Event2NLL

# Convert catalogue to NLL format
nll = Event2NLL(catalog=catalog,  # Previously derived seismicity catalog
                nll_basepath="../NonLinLoc",
                vel_model=vel_model,  # Velocity model for NLL
                stations_df=stations, # pandas Dataframe for stations
                )
```

The `.obs` files are saved in `nll_basepath/obs` (e.g. `../NonLinLoc/obs/`).
From here, the files can be copied to a different NonLinLoc project for relocalisation.

### Notes on NonLinLoc Implementation
Note the NLL implementation can be modified in `association.core.nll_functions` since
some settings for NLL are hard coded.

## Examples
- `tests/example_harpa.py`: Association using HARPA.
- `tests/example_pyocto.py`: Association using PyOcto (including velocity model).
- `tests/example_gamma.py`: Association using GaMMA
- `tests/example_harpa_simulations.py`: Simulates events in area, adds random picks and associates
   simulated picks with HARPA.
- `tests/compare_associators.py`: Comparison of all associators classes against each other on
   example dataset.
