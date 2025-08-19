import json

import obspy
import pandas as pd
import numpy as np
from typing import Union

from pyproj import Proj
from concurrent.futures import ThreadPoolExecutor
from harpa import association  # noqa
from obspy.core.event.origin import Pick, Origin
from obspy.core.event.event import Event
from obspy.core.event.base import WaveformStreamID


def convert_station_json(stations: dict) -> pd.DataFrame:
    """

    :param stations:
    :return:
    """
    station_df = []
    for station in stations.keys():
        station_df.append(
            {
                "id": station,
                "latitude": stations[station]["coords"][0],
                "longitude": stations[station]["coords"][1],
                "elevation": stations[station]["coords"][2],
            }
        )

    return pd.DataFrame(station_df)


def load_stations(station_json: str):
    """

    :param station_json:
    :return:
    """
    with open(station_json) as f_json:
        stations = json.load(f_json)
        stations = convert_station_json(stations)

    return stations


def area_limits(stations: (pd.DataFrame, dict)) -> dict:
    """

    :param stations:
    :return:
    """
    # Convert stations to pandas dataframe
    if isinstance(stations, dict):
        stations = convert_station_json(stations=stations)

    limits = {
        "latitude": (
            np.min(stations["latitude"] - 0.2),
            np.max(stations["latitude"]) + 0.2,
        ),
        "longitude": (
            np.min(stations["longitude"] - 0.2),
            np.max(stations["longitude"]) + 0.2,
        ),
    }

    # Add center to limits
    limits["center"] = (
        (limits["latitude"][1] - limits["latitude"][0]) / 2 + limits["latitude"][0],
        (limits["longitude"][1] - limits["longitude"][0]) / 2 + limits["longitude"][0],
    )

    return limits


def sort_events(events: list):
    # Create list with all dates
    dates = [event.origins[0].time.datetime for event in events]

    zipped_pairs = zip(dates, events)
    try:
        sorted_events = [x for _, x in sorted(zipped_pairs)]
    except TypeError:
        print(zipped_pairs)
        print("# BEGIN MESSAGE #\nDid not sort events\n# END MESSAGE #")
        return events

    return sorted_events


def associate_harpa(
    picks: pd.DataFrame,
    config: dict,
    station_json: Union[str, pd.DataFrame],
    verbose: int = 10,
    iterations: int = 3,
    second_pass: bool = False,
) -> obspy.Catalog:
    if isinstance(station_json, str) is True:
        station_df = load_stations(station_json=station_json)
    else:
        station_df = station_json

    x0 = station_df["longitude"].median()
    y0 = station_df["latitude"].median()
    xmin = station_df["longitude"].min()
    xmax = station_df["longitude"].max()
    ymin = station_df["latitude"].min()
    ymax = station_df["latitude"].max()

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
    station_df[["x(km)", "y(km)"]] = station_df.apply(
        lambda x: pd.Series(proj(longitude=x.longitude, latitude=x.latitude)), axis=1
    )
    station_df["z(km)"] = -station_df["elevation"] / 1e3

    # Start association with HARPA
    with ThreadPoolExecutor(max_workers=8) as executor:
        future = executor.submit(association, picks, station_df, config, verbose)
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

    event_lst = []
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
        second_pass_events = associate_harpa(
            picks=remaining_picks,
            config=config,
            station_json=station_df,
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


def double(x, iterations=2):
    for i in range(iterations):
        x = double(x, iterations=0)
    return 2 * x


if __name__ == "__main__":
    print(double(2))
