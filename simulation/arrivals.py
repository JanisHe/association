"""
Simulate arrival times for stations and random locations within the seismic network
TODO: Allow single phase arrivals in create event
"""

import random
import pandas as pd
import numpy as np
from typing import Optional

from obspy.geodetics import degrees2kilometers, locations2degrees
from obspy import UTCDateTime


def create_event(
    stations: pd.DataFrame,
    latitude_range: tuple,
    longitude_range: tuple,
    depth_range: tuple,
    origin_time: UTCDateTime,
    vp: float = 4500,
    vs: Optional[float] = None,
    min_num_stations: int = 4,
) -> (pd.DataFrame, tuple):
    if not vs:
        vs = vp / np.sqrt(3)

    # Simulate location of random earthquake
    latitude = np.random.uniform(low=min(latitude_range), high=max(latitude_range))
    longitude = np.random.uniform(low=min(longitude_range), high=max(longitude_range))
    depth = np.random.uniform(low=min(depth_range), high=max(depth_range))

    # Select stations that have a pick by choosing station indices from all stations
    station_idx = []
    num_stations = np.random.randint(min_num_stations, len(stations))
    all_idx = list(stations.index)
    for i in range(num_stations):
        idx = np.random.choice(all_idx)
        all_idx.remove(idx)
        station_idx.append(idx)

    # Create arrivals for each station
    arrivals = {"trace_id": [], "peak_time": [], "phase": [], "peak_value": []}
    for stat_idx in station_idx:
        # Compute traveltimes
        distance_degree = locations2degrees(
            lat1=stations.loc[stat_idx, "latitude"],
            long1=stations.loc[stat_idx, "longitude"],
            lat2=latitude,
            long2=longitude,
        )
        epi_distance_m = degrees2kilometers(degrees=distance_degree) * 1000
        hypo_distance_m = np.sqrt(
            (stations.loc[stat_idx, "elevation_m"] + depth) ** 2 + epi_distance_m**2
        )

        for phase, v in zip(["P", "S"], [vp, vs]):
            traveltime = hypo_distance_m / v
            arrivals["phase"].append(phase)
            arrivals["trace_id"].append(stations.loc[stat_idx, "trace_id"])
            arrivals["peak_time"].append(origin_time + traveltime)
            arrivals["peak_value"].append(np.random.uniform(low=0.2, high=1.0))

    return pd.DataFrame(arrivals), (latitude, longitude, depth)


def add_random_picks(event_df: pd.DataFrame, percentage: float = 15) -> pd.DataFrame:
    """
    Adding fake picks to event_df

    :param event_df:
    :param percentage:
    :return:
    """
    stations = list(
        set(event_df["trace_id"])
    )  # Create list with all stations from event_df
    starttime, endtime = (
        min(event_df["peak_time"]) - 30,
        max(event_df["peak_time"]) + 30,
    )  # Select time range from all times (add 30 seconds)
    add_picks = int(
        len(event_df) * percentage / 100
    )  # Number of how many events will be added

    # Transform pandas dataframe to dictionary
    event_dct = event_df.to_dict(orient="list")

    # Create random picks with low probabilities
    for i in range(add_picks):
        event_dct["trace_id"].append(random.choice(stations))
        event_dct["peak_time"].append(
            starttime + np.random.uniform(low=0, high=endtime - starttime)
        )
        event_dct["phase"].append(random.choice(["P", "S"]))
        event_dct["peak_value"].append(np.random.uniform(low=0.2, high=0.5))

    return pd.DataFrame(event_dct)


if __name__ == "__main__":
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

    for i in range(10):
        origin_time = starttime + np.random.randint(0, 300)
        df, _ = create_event(
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

    df_events = add_random_picks(event_df=df_events)
    print()
