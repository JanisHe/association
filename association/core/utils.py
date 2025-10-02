import pandas as pd


def area_limits(stations: pd.DataFrame, lat_lon_eps: float = 0.2) -> dict:
    limits = {
        "latitude": (
            min(stations["latitude"]) - lat_lon_eps,
            max(stations["latitude"]) + lat_lon_eps,
        ),
        "longitude": (
            min(stations["longitude"]) - lat_lon_eps,
            max(stations["longitude"]) + lat_lon_eps,
        ),
    }

    # Compute center of area and add to limits
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
