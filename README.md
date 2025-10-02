# Seismic Phase Association for AIS Project

This packages contains scripts and interfaces to run either [Pyocto](https://github.com/yetinam/pyocto)
or [Harpa](https://github.com/DaDaCheng/phase_association/tree/main) for seismic phase association on PhaseNet picks.

## Required packages
* standard packages: `numpy`, `obspy`, `pandas`, `seisbench`, `matplotlib`
* Packages for phase association:
  - `pyocto`: `pip install pyocto`,
     PyOcto requires Pyrocko: `https://github.com/pyrocko/pyrocko`
  - `harpa`: `pip install -q git+https://github.com/DaDaCheng/phase_association.git`,
     Harpa requires POT: `pip install POT`
     Harpa requires scikit-learn: `pip install scikit-learn`

## Preparing files
### Stations
A `.csv` file that contains station information (note the file header!):

| id         | latitude | longitude | elevation_m |
|------------|----------|-----------|-------------|
| FO.BETS.00 | 48.89357 | 7.92429   | 146         |
| RG.RITT.00 | 48.89436 | 7.96103   | 138         |
| RG.KUHL.00 | 48.91473 | 7.92996   | 176         |
| FO.OPS.00  | 48.92126 | 7.88278   | 198         |

- `id`: ID of each station (`network.station.location`)
- `latitude`: Latitude of station in degree
- `longitude`: Longitude of station in degree
- `elevation_m`: Elevation of station in m above sea level

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

- `id`: ID of station (`network.station.location`)
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
Once all files are created, phase association is done by running the following code
(config dictionaries can be varied by user; use available keywords that can be found in
the descriptions of [Pyocto](https://github.com/yetinam/pyocto)
or [Harpa](https://github.com/DaDaCheng/phase_association/tree/main)):
```
from association.core.interfaces import interface_pyocto, interface_harpa

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

Both interface functions return an `obspy.Catalog` that contains all events, including
the picks.

## Relocalisation with NonLinLoc
After seismic phase association, the events can be relocalised using [NonLinLoc](http://alomax.free.fr/nlloc/).
Therefore, a velocity model (see above) is required.

```
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

## Examples
- `tests/example_harpa.py`: Association using HARPA.
- `tests/example_pyocto.py`: Association using PyOcto.
- `tests/example_harpa_simulations.py`: Simulates events in area, adds random picks and associates
   simulated picks with HARPA.
