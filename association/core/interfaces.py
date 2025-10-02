"""
Function for different association algorithms:
HARPA, PyOcto, and GaMMA
"""

import os
import datetime
import pyocto
import obspy
import copy

import pandas as pd
from obspy import UTCDateTime
from pyproj import Proj

from typing import Optional
from concurrent.futures import ThreadPoolExecutor
from obspy.geodetics import degrees2kilometers
from obspy.core.event.origin import Pick, Origin
from obspy.core.event.event import Event
from obspy.core.event.base import WaveformStreamID

from harpa import association  # noqa
from association.core.utils import area_limits, sort_events
from gamma.utils import association as association_gamma


def interface_gamma(
    picks: pd.DataFrame,
    config: dict,
    stations: pd.DataFrame,
    velocity_model: Optional[pd.DataFrame] = None,
    verbose: bool = False,
) -> obspy.Catalog:
    """

    :param picks:
    :param config:
    :param stations:
    :param velocity_model:
    :param verbose:
    :return:
    """
    # Get grid from latitude and longitude
    limits = area_limits(stations=stations)
    config["xlim_degree"] = limits["latitude"]
    config["ylim_degree"] = limits["longitude"]

    # Create projection from pyproj. Use a string that is then transfered to Proj (https://proj.org)
    proj = Proj(
        f"+proj=sterea +lat_0={limits['center'][0]} +lon_0={limits['center'][1]} +units=km"
    )

    # Compute boundaries for x and y
    ylim1, xlim1 = proj(
        latitude=config["xlim_degree"][0], longitude=config["ylim_degree"][0]
    )
    ylim2, xlim2 = proj(
        latitude=config["xlim_degree"][1], longitude=config["ylim_degree"][1]
    )
    config["x(km)"] = [xlim1, xlim2]
    config["y(km)"] = [ylim1, ylim2]

    # Make config dict for GaMMA
    # Config for earthquake location
    config["ncpu"] = config.get("ncpu") or os.cpu_count()
    config["dims"] = ["x(km)", "y(km)", "z(km)"]
    config["use_dbscan"] = config.get("use_dbscan") or True  # GaMMA default is True
    config["use_amplitude"] = (
        config.get("use_amplitude") or False
    )  # Set to False if not chosen by user
    config["z(km)"] = config.get("z(km)")
    config["method"] = config.get("method") or "BGMM"  # GMM is also possible
    config["initial_points"] = config.get("inital_points") or [
        1,
        1,
        1,
    ]  # Use GaMMA default value
    config["covariance_prior"] = config.get("covariance_prior") or [
        5,
        5,
    ]  # Use GaMMA default value

    # Build velocity model
    if isinstance(velocity_model, pd.DataFrame):
        config["eikonal"] = {
            "vel": {
                "p": velocity_model["vp"].to_list(),
                "s": velocity_model["vs"].to_list(),
                "z": velocity_model["depth"].to_list(),
            },
            "h": 1.0,
            "xlim": config["x(km)"],
            "ylim": config["y(km)"],
            "zlim": config["z(km)"],
        }
        if config.get("vel"):
            config.pop("vel")  # Remove vel from config if velocity model is given
    else:
        config["vel"] = config.get("vel") or {
            "p": 6.0,
            "s": 6.0 / 1.75,
        }  # Use default GaMMA value

    config["bfgs_bounds"] = (
        (config["x(km)"][0] - 1, config["x(km)"][1] + 1),
        (config["y(km)"][0] - 1, config["y(km)"][1] + 1),
        (0, config["z(km)"][1] + 1),  # x
        (None, None),  # t
    )

    # The initial number of clusters is determined by
    # (Number of picks)/(Number of stations) * (oversampling factor).
    config["oversample_factor"] = config.get("oversample_factor") or 10  # Default is 10

    # DBSCAN
    config["dbscan_eps"] = config.get("dbscan_eps") or 10  # Default is 10 s
    config["dbscan_min_samples"] = config.get("dbscan_min_samples") or 3  # Default is 3

    # Filtering
    config["min_picks_per_eq"] = config.get("min_picks_per_eq")
    config["min_p_picks_per_eq"] = config.get("min_p_picks_per_eq")
    config["min_s_picks_per_eq"] = config.get("min_s_picks_per_eq")
    config["max_sigma11"] = config.get("max_sigma11")
    config["max_sigma22"] = config.get("max_sigma22") or 1.0
    config["max_sigma12"] = config.get("max_sigma12")

    # Print configs if verbose is True
    if verbose:
        print("Config for GaMMA:")
        for key, value in config.items():
            print(f"{key}: {value}")
        print()

    # Transform station coordinates
    stations[["x(km)", "y(km)"]] = stations.apply(
        lambda x: pd.Series(proj(longitude=x.longitude, latitude=x.latitude)), axis=1
    )
    stations["z(km)"] = stations["elevation"].apply(lambda depth: -depth / 1e3)

    # Do association
    catalogs, assignments = association_gamma(
        picks, stations, config, method=config["method"]
    )

    catalog = pd.DataFrame(catalogs)
    assignments = pd.DataFrame(
        assignments, columns=["pick_index", "event_index", "gamma_score"]
    )
    if len(catalog) > 0:
        catalog.sort_values(by=["time"], inplace=True, ignore_index=True)
    else:
        return obspy.Catalog([])

    # Transform earthquake locations to lat long
    event_lat = []
    event_long = []
    for x, y in zip(catalog["x(km)"], catalog["y(km)"]):
        long, lat = proj(x, y, inverse=True)
        event_lat.append(lat)
        event_long.append(long)
    catalog["latitude"] = event_lat
    catalog["longitude"] = event_long

    # Create earthquake catalog
    event_lst = []  # Empty event list to save obspy events
    # Loop over each event in GaMMA catalog
    for event_count in range(len(catalog)):
        event_idx = catalog["event_index"][event_count]
        event_picks = [
            picks.iloc[i]
            for i in assignments[assignments["event_index"] == event_idx]["pick_index"]
        ]

        # Create dictionary for each station that contains P and S phases
        picks_lst = []  # Create list with obspy picks class
        for p in event_picks:
            waveform_id = WaveformStreamID(
                network_code=p["id"].split(".")[0],
                station_code=p["id"].split(".")[1],
                location_code=p["id"].split(".")[2],
            )
            pick = Pick(
                time=obspy.UTCDateTime(p["timestamp"]),
                waveform_id=waveform_id,
                phase_hint=p["type"].upper(),
            )
            picks_lst.append(pick)

        origin = Origin(
            time=catalog["time"][event_count],
            longitude=catalog["longitude"][event_count],
            latitude=catalog["latitude"][event_count],
            depth=catalog["z(km)"][event_count] * 1e3,
        )

        event = Event(
            picks=picks_lst,
            force_resource_id=True,
            origins=[origin],
        )

        # Append event to final event list
        event_lst.append(event)

    return obspy.Catalog(event_lst)


def interface_harpa(
    picks: pd.DataFrame,
    config: dict,
    stations: pd.DataFrame,
    verbose: int = 10,
    iterations: int = 3,
    second_pass: bool = False,
) -> obspy.Catalog:
    """

    :param picks:
    :param config:
    :param stations:
    :param verbose:
    :param iterations:
    :param second_pass:
    :return:
    """

    x0 = stations["longitude"].mean()
    y0 = stations["latitude"].mean()
    xmin = stations["longitude"].min()
    xmax = stations["longitude"].max()
    ymin = stations["latitude"].min()
    ymax = stations["latitude"].max()

    # Update config
    config["center"] = (x0, y0)
    config["xlim_degree"] = (2 * xmin - x0, 2 * xmax - x0)
    config["ylim_degree"] = (2 * ymin - y0, 2 * ymax - y0)

    proj = Proj(
        f"+proj=sterea +lon_0={config['center'][0]} +lat_0={config['center'][1]} +units=km"
    )
    config["x(km)"] = proj(
        longitude=config["xlim_degree"], latitude=[config["center"][1]] * 2
    )[0]
    config["y(km)"] = proj(
        longitude=[config["center"][0]] * 2, latitude=config["ylim_degree"]
    )[1]
    stations[["x(km)", "y(km)"]] = stations.apply(
        lambda x: pd.Series(proj(longitude=x.longitude, latitude=x.latitude)), axis=1
    )
    stations["z(km)"] = -stations["elevation"] / 1e3

    # Start association with HARPA
    with ThreadPoolExecutor(max_workers=8) as executor:
        future = executor.submit(association, picks, stations, config, verbose)
    pick_df_out, catalog_df = future.result()

    # Create obspy catalogue
    # Transform earthquake locations to lat long
    event_lat = []
    event_long = []
    for x, y in zip(catalog_df["x(km)"], catalog_df["y(km)"]):
        long, lat = proj(x, y, inverse=True)
        event_lat.append(lat)
        event_long.append(long)
    catalog_df["latitude"] = event_lat
    catalog_df["longitude"] = event_long

    event_lst = []  # Initialize empty list for events
    # Loop over each event in HARPA catalog
    for event_count in range(len(catalog_df)):
        event_idx = catalog_df["event_index"][event_count]
        event_picks_df = pick_df_out[pick_df_out["event_index"] == event_idx]
        # event_picks = [pick_df_out.iloc[i] for i in range(pick_df_out[pick_df_out["event_index"] == event_idx])]

        # Create dictionary for each station that contains P and S phases
        # Create list with obspy picks class
        picks_lst = []
        for idx in list(event_picks_df.index):
            waveform_id = WaveformStreamID(
                network_code=event_picks_df.loc[idx, "station_id"].split(".")[0],
                station_code=event_picks_df.loc[idx, "station_id"].split(".")[1],
                location_code=event_picks_df.loc[idx, "station_id"].split(".")[2],
            )
            pick = Pick(
                time=obspy.UTCDateTime(event_picks_df.loc[idx, "timestamp"]),
                waveform_id=waveform_id,
                phase_hint=event_picks_df.loc[idx, "type"].upper(),
            )
            picks_lst.append(pick)

        origin = Origin(
            time=catalog_df["time"][event_count],
            longitude=catalog_df["longitude"][event_count],
            latitude=catalog_df["latitude"][event_count],
            depth=catalog_df["z(km)"][event_count] * 1e3,
        )
        ev = Event(picks=picks_lst, force_resource_id=True, origins=[origin])

        # Append event to final event list
        event_lst.append(ev)

    # Run HARPA a second time on no used picks
    # Picks that have not been used are marked with the event_index = -1 in pick_df_out
    if iterations > 1:
        remaining_picks = pick_df_out[pick_df_out["event_index"] == -1]
        # Modify remaining picks dataframe
        remaining_picks = remaining_picks.rename(columns={"station_id": "id"})
        second_pass_events = interface_harpa(
            picks=remaining_picks,
            config=config,
            stations=stations,
            verbose=verbose,
            second_pass=True,
            iterations=iterations - 1,
        )
        event_lst += second_pass_events

    if second_pass is True:
        return event_lst
    else:
        # Sort event_lst
        event_lst = sort_events(events=event_lst)
        return obspy.Catalog(events=event_lst)


def interface_pyocto(
    stations: pd.DataFrame,
    picks: pd.DataFrame,
    config: dict,
    velocity_model: Optional[pd.DataFrame] = None,
    station_column_renaming: Optional[dict[str, str]] = None,
    pick_column_renaming: Optional[dict[str, str]] = None,
) -> obspy.Catalog:
    """

    :param stations:
    :param picks:
    :param config
    :param velocity_model:
    :param station_column_renaming: {"trace_id": "id", "elevation_m": "elevation"}
    :param pick_column_renaming: {"id": "station", "timestamp": "time", "type": "phase"}
    :return:

    Required parameters: zlim: tuple[int, int]
    Required if velocity_model is used: velocity_model_path: pathname to save velocity model
    Optional if velocity_model is used: tolerance, delta

    Required if no velocity_model is used: p_velocity: constant P-velocity for model
                                           s_velocity: constant S-velocity for model
    Optional if no velocity_model is used: tolerance

    Other parameters/keyword args will be handed over to the PyOcto associator class and will be checked there.
    """
    config = copy.deepcopy(config)  # Create copy of conifg to avoid overwriting
    area = area_limits(stations=stations)  # Get limits and center of area

    # Check if velocity_model_path is in config if velocity model is used
    if isinstance(velocity_model, pd.DataFrame):
        try:
            config["velocity_model_path"]
        except KeyError:
            msg = (
                "Key 'velocity_model_path' is missing in config to save velocity model."
            )
            raise ValueError(msg)

    # Build velocity model (if available)
    max_degree = max(
        [
            abs(max(area["latitude"]) - min(area["latitude"])),
            abs(max(area["longitude"]) - min(area["longitude"])),
        ]
    )
    if isinstance(velocity_model, pd.DataFrame):
        pyocto.VelocityModel1D.create_model(
            model=velocity_model,
            delta=config["delta"] if config.get("delta") else 1,
            xdist=degrees2kilometers(degrees=max_degree),
            zdist=max(config["zlim"]),
            path=config["velocity_model_path"],
        )
        velocity_model = pyocto.VelocityModel1D(
            path=config.pop("velocity_model_path"),
            tolerance=config["tolerance"] if config.get("tolerance") else 2.0,
        )
    else:  # Constant velocity model
        velocity_model = pyocto.VelocityModel0D(
            p_velocity=config.pop("p_velocity"),
            s_velocity=config.pop("s_velocity"),
            tolerance=config["tolerance"] if config.get("tolerance") else 2.0,
        )

    # Set up associator from PyOcto
    associator = pyocto.OctoAssociator.from_area(
        lat=(min(stations["latitude"]), max(stations["latitude"])),
        lon=(min(stations["longitude"]), max(stations["longitude"])),
        velocity_model=velocity_model,
        **config,
    )

    # Convert stations to PyOcto format (i.e. required column names are 'id', 'latitude',
    # 'longitude', 'elevation')
    if station_column_renaming:
        stations = stations.rename(columns=station_column_renaming)

    # Convert picks to PyOcto format (i.e. required column names are 'time', 'station', 'phase')
    if pick_column_renaming:
        picks = picks.rename(columns=pick_column_renaming)

    picks["time"] = picks["time"].apply(
        lambda x: UTCDateTime(x)
    )  # Convert time strings to obspy UTCDateTime
    picks["time"] = picks["time"].apply(
        lambda x: x.datetime.timestamp()
    )  # Reformat time
    picks["phase"] = picks["phase"].apply(
        lambda x: x.upper()
    )  # Reformat phase to upper case letters

    # Convert stations to required format
    stations = associator.transform_stations(stations=stations)

    # Start association
    events, assignments = associator.associate(picks, stations)
    associator.transform_events(events)

    # Proof whether events have been detected
    if len(events) == 0:
        return obspy.Catalog()

    # Assign correct time to event
    events["time"] = events["time"].apply(datetime.datetime.fromtimestamp)
    assignments["time"] = assignments["time"].apply(datetime.datetime.fromtimestamp)

    # Merge events to pick assignments
    assignments = pd.merge(
        events, assignments, left_on="idx", right_on="event_idx", suffixes=("", "_pick")
    )

    # Create final obspy catalogue from events and picks
    # Picks are stored in a dictionary for each event by idx of event
    picks_dct = {}
    origins = {}
    for i in range(len(assignments)):
        event_idx = assignments.loc[i, "idx"]
        if event_idx not in picks_dct.keys():
            picks_dct.update({event_idx: []})  # List with all picks for an event

        if event_idx not in origins.keys():
            origins.update(
                {
                    event_idx: Origin(
                        time=obspy.UTCDateTime(assignments.loc[i, "time"]),
                        latitude=assignments.loc[i, "latitude"],
                        longitude=assignments.loc[i, "longitude"],
                        depth=assignments.loc[i, "depth"] * 1000,
                    )
                }
            )
        # Create obspy pick
        network_code, station_code, location = assignments.loc[i, "station"].split(".")
        waveform_id = WaveformStreamID(
            network_code=network_code, station_code=station_code, location_code=location
        )
        pick = Pick(
            time=obspy.UTCDateTime(assignments.loc[i, "time_pick"]),
            waveform_id=waveform_id,
            phase_hint=assignments.loc[i, "phase"],
        )
        pick.time_errors["uncertainty"] = abs(assignments.loc[i, "residual"])

        # Append pick to event_idx in dictionary
        picks_dct[event_idx].append(pick)

    # Create obspy catalogue from picks_dct and origins
    event_list = []
    for event_idx in origins.keys():
        event_list.append(
            Event(
                origins=[origins[event_idx]],
                picks=picks_dct[event_idx],
                force_resource_id=True,
            )
        )

    # Sort events by date
    event_list = sort_events(events=event_list)

    return obspy.Catalog(event_list)


def association_wrapper(
    stations: pd.DataFrame, picks: pd.DataFrame, method: str, **kwargs
):
    pass
