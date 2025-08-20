"""
Function for different association algorithms:
HARPA, PyOcto, and GaMMA
"""
import datetime
import pyocto
import obspy

import pandas as pd

from typing import Optional
from obspy.geodetics import degrees2kilometers
from obspy.core.event.origin import Pick, Origin
from obspy.core.event.event import Event
from obspy.core.event.base import WaveformStreamID

from core.utils import area_limits, sort_events


# TODO: Input: csv file with picks, output catalog with events
#       Columns of csv file: trace_id (network.station.location), peak_time (UTCDateTime ???), peak_value, phase
#       csv file needs to be created before (i.e., concatenating files from all stations and consider overlapping in time)


def interface_pyocto(stations: pd.DataFrame,
                     picks: pd.DataFrame,
                     velocity_model: Optional[pd.DataFrame] = None,
                     **pyocto_kwargs) -> obspy.Catalog:
    """

    :param stations:
    :param picks:
    :param velocity_model:
    :param pyocto_kwargs:
    :return:

    List of keywords: delta, zdist, velocity_model_path
    """

    area = area_limits(stations=stations)  # Get limits and center of area

    # Build velocity model (if available)
    max_degree = max([abs(max(area["latitude"]) - min(area["latitude"])),
                      abs(max(area["longitude"]) - min(area["longitude"]))])
    if velocity_model:
        pyocto.VelocityModel1D.create_model(model=velocity_model,
                                            delta=pyocto_kwargs["delta"] if pyocto_kwargs.get("delta") else 1,
                                            xdist=degrees2kilometers(degrees=max_degree),
                                            zdist=max(pyocto_kwargs["zdist"]),
                                            path=pyocto_kwargs["velocity_model_path"])
        velocity_model = pyocto.VelocityModel1D(path=pyocto_kwargs["velocity_model_path"],
                                                tolerance=pyocto_kwargs["tolerance"] if pyocto_kwargs.get("tolerance") else 2.0)
    else:  # Constant velocity model
        velocity_model = pyocto.VelocityModel0D(
            p_velocity=7.0,
            s_velocity=4.0,
            tolerance=2.0,
        )


    # Set up associator from PyOcto
    associator = pyocto.OctoAssociator.from_area(
        lat=(min(stations["latitude"]), max(stations["latitude"])),
        long=(min(stations["longitude"]), max(stations["longitude"])),
        velocity_model=velocity_model,
        **pyocto_kwargs
    )

    # Convert stations to required format
    stations = associator.transform_stations(stations=stations)

    # Start association
    events, assignments = associator.associate_seisbench(picks, stations)
    associator.transform_events(events)

    # Proof whether events have been detected
    if len(events) == 0:
        return obspy.Catalog()

    # Assign correct time to event
    events["time"] = events["time"].apply(datetime.datetime.fromtimestamp, tz=datetime.timezone.utc)
    assignments["time"] = assignments["time"].apply(datetime.datetime.fromtimestamp, tz=datetime.timezone.utc)

    # Merge events to pick assignments
    assignments = pd.merge(events, assignments, left_on="idx", right_on="event_idx", suffixes=("", "_pick"))

    # Create final obspy catalogue from events and picks
    # Picks are stored in a dictionary for each event by idx of event
    picks_dct = {}
    origins = {}
    for i in range(len(assignments)):
        event_idx = assignments.loc[i, "idx"]
        if event_idx not in picks_dct.keys():
            picks_dct.update({event_idx: []})  # List with all picks for an event

        if event_idx not in origins.keys():
            origins.update({event_idx: Origin(time=obspy.UTCDateTime(assignments.loc[i, "time"]),
                                              latitude=assignments.loc[i, "latitude"],
                                              longitude=assignments.loc[i, "longitude"],
                                              depth=assignments.loc[i, "depth"] * 1000)
                            })
        # Create obspy pick
        network_code, station_code, location = assignments.loc[i, "station"].split(".")
        waveform_id = WaveformStreamID(network_code=network_code,
                                       station_code=station_code,
                                       location_code=location)
        pick = Pick(time=obspy.UTCDateTime(assignments.loc[i, "time_pick"]),
                    waveform_id=waveform_id,
                    phase_hint=assignments.loc[i, "phase"]
                    )
        pick.time_errors["uncertainty"] = abs(assignments.loc[i, "residual"])

        # Append pick to event_idx in dictionary
        picks_dct[event_idx].append(pick)

    # Create obspy catalogue from picks_dct and origins
    event_list = []
    for event_idx in origins.keys():
        event_list.append(Event(origins=[origins[event_idx]],
                                picks=picks_dct[event_idx],
                                force_resource_id=True)
                          )

    # Sort events by date
    event_list = sort_events(events=event_list)

    return obspy.Catalog(event_list)


def interface_gamma():
    pass


def interface_harpa():
    pass

def association_wrapper(stations: pd.DataFrame,
                        picks: pd.DataFrame,
                        method: str,
                        **kwargs):
    pass
