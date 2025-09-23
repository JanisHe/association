from typing import Optional

import pathlib
import obspy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def moveout(
        catalog: obspy.Catalog,
        stations: pd.DataFrame,
        ax_lat: Optional[plt.axes] = None,
        ax_lon: Optional[plt.axes] = None,
        csv_pick_files: Optional[list] = None,
        title: str = ""
):
    """
    Plotting moveout curve from catalog along time and latitude / longitude
    """
    # Set color for each line
    colors = [
        "tab:blue",
        "tab:orange",
        "tab:green",
        "tab:red",
        "tab:purple",
        "tab:brown",
        "tab:pink",
        "tab:gray",
        "tab:olive",
        "tab:cyan",
    ]

    # Set title for ax_lat and ax_lon
    if ax_lon and len(title) > 0:
        ax_lon.set_title(title)
    if ax_lat and len(title) > 0:
        ax_lat.set_title(title)

    # Estimate number of needed colors
    colors = colors * int(np.ceil(len(catalog) / len(colors)))

    # If pick_files available plot all available picks in grey
    csv_p_picks = []
    csv_s_picks = []
    csv_p_latitudes = []
    csv_s_latitudes = []
    csv_p_longitudes = []
    csv_s_longitudes = []
    if csv_pick_files:
        for filename in csv_pick_files:
            # Find location of station from station_json
            station_code = pathlib.Path(filename).stem
            try:
                station_idx = list(stations["id"]).index(station_code)
            except ValueError:
                continue
            picks = pd.read_csv(filename)  # Read pick-file

            # Loop over each entry in csv-file
            for idx in range(len(picks)):
                if picks.loc[idx, "phase"].lower() == "p":
                    csv_p_picks.append(
                        obspy.UTCDateTime(picks.loc[idx, "peak_time"]).datetime
                    )
                    csv_p_latitudes.append(stations.loc[station_idx, "latitude"])
                    csv_p_longitudes.append(stations.loc[station_idx, "longitude"])
                elif picks.loc[idx, "phase"].lower() == "s":
                    csv_s_picks.append(
                        obspy.UTCDateTime(picks.loc[idx, "peak_time"]).datetime
                    )
                    csv_s_latitudes.append(stations.loc[station_idx, "latitude"])
                    csv_s_longitudes.append(stations.loc[station_idx, "longitude"])

        # Plot picks in grey
        if ax_lat:
            ax_lat.scatter(
                x=csv_p_picks,
                y=csv_p_latitudes,
                color="grey",
                alpha=0.35,
                rasterized=True,
            )
        if ax_lon:
            ax_lon.scatter(
                x=csv_p_picks,
                y=csv_p_longitudes,
                color="grey",
                alpha=0.35,
                rasterized=True,
            )

        # Plot S-picks
        if ax_lat:
            ax_lat.scatter(
                x=csv_s_picks,
                y=csv_s_latitudes,
                color="grey",
                alpha=0.35,
                edgecolors=["k"] * len(csv_s_picks),
                rasterized=True,
            )
        if ax_lon:
            ax_lon.scatter(
                x=csv_s_picks,
                y=csv_s_longitudes,
                color="grey",
                alpha=0.35,
                edgecolors=["k"] * len(csv_s_picks),
                rasterized=True,
            )

    for idx, event in enumerate(catalog):
        # origin_time = event.origins[-1].datetime
        p_picks = []
        s_picks = []
        p_latitudes = []
        p_longitudes = []
        s_latitudes = []
        s_longitudes = []
        for pick in event.picks:
            # Find latitude and longitude of station
            station_code = (
                f"{pick.waveform_id.network_code}.{pick.waveform_id.station_code}."
                f"{pick.waveform_id.location_code}"
            )
            try:
                station_idx = list(stations["id"]).index(station_code)
            except ValueError:
                continue

            if pick.phase_hint.lower() in ["p", "pg", "pn"]:
                p_picks.append(pick.time.datetime)
                p_latitudes.append(stations.loc[station_idx, "latitude"])
                p_longitudes.append(stations.loc[station_idx, "longitude"])
            elif pick.phase_hint.lower() in ["s", "sg", "sn"]:
                s_picks.append(pick.time.datetime)
                s_latitudes.append(stations.loc[station_idx, "latitude"])
                s_longitudes.append(stations.loc[station_idx, "longitude"])

        # Sort latitudes
        # zipped = zip(p_latitudes, p_picks)

        # Plot moveout
        # ax.plot(p_picks,
        #         p_latitudes,
        #         color="grey",
        #         linestyle="--",
        #         linewidth=0.5)

        # Plot P-picks
        if ax_lat:
            ax_lat.scatter(x=p_picks, y=p_latitudes, color=colors[idx], rasterized=True)
        if ax_lon:
            ax_lon.scatter(
                x=p_picks, y=p_longitudes, color=colors[idx], rasterized=True
            )

        # Plot S-picks
        if ax_lat:
            ax_lat.scatter(
                x=s_picks,
                y=s_latitudes,
                color=colors[idx],
                edgecolors=["k"] * len(s_picks),
                rasterized=True,
            )
        if ax_lon:
            ax_lon.scatter(
                x=s_picks,
                y=s_longitudes,
                color=colors[idx],
                edgecolors=["k"] * len(s_picks),
                rasterized=True,
            )
