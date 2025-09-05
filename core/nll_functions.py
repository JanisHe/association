import os
import glob
import shutil
import datetime
import obspy

import pandas as pd
from typing import Union, Optional

from obspy.core.event import read_events, Catalog

from harpa.helpers import load_stations, area_limits
from core.utils import sort_events


class Event2NLL:
    """
    Class to relocate events with NonLinLoc.
    The class takes an obspy catalog as input, which was previously created using a phase association algorithm.
    The catalog must contain the seismic phase onset. Afterwards, the relocation is done with NonLinLoc.

    :param catalog: Catalog to relocalise data
    :param nll_basepath: Full pathname to save results from NonLinLoc
    :param vel_model: Full pathname of velocity model
    :param stations_df: csv-filename or pandas dataframe that contains station information
    :param nll_executable: Executable of NonLinLoc
    :param create_files: Set to True if all required NLL files should be created. Default is True.
                         Note, existing files will be overwritten
    :param pick_uncertainty: Assign uncertainty to picks, if not available from picks in catalog. Default is 0.5 s.
    :param kwargs: Possible kwargs: 'signature' and 'comment'.
    """
    def __init__(self,
                 catalog: Catalog,
                 nll_basepath: str,
                 vel_model: str,
                 stations_df: Union[str, pd.DataFrame],
                 nll_executable: str = "/work/software/nlloc/7.00.16/src/bin",
                 create_files: bool = True,
                 pick_uncertainty: float = 0.5,
                 **kwargs):
        self.catalog = catalog
        self.nll_basepath = nll_basepath
        self.vel_model = vel_model
        self.nll_executables = nll_executable
        self.stations_json = stations_df
        self.pick_uncertainty = pick_uncertainty

        self.event_names = [str(event.origins[0].time) for event in catalog]
        self.config_name = "_".join(str(datetime.datetime.now()).split())
        self.df_vel_model = None
        self.conf_file = os.path.join(self.nll_basepath, 'conf', f'{self.config_name}.conf')

        # Read stations
        if isinstance(stations_df, pd.DataFrame):
            self.df_stations = stations_df
        else:
            self.df_stations = load_stations(station_json=stations_df)

        # Compute midpoint from all stations in station_json
        limits = area_limits(stations=self.df_stations)
        self.latorig = limits["center"][0]
        self.longorig = limits["center"][1]

        if create_files is True:
            self.create_nll_dirs()
            self.create_obs()
            self.create_config(**kwargs)
        else:
            self.conf_file = glob.glob(os.path.join(self.nll_basepath, 'conf', '*'))[0]

    def create_nll_dirs(self,
                        delete_basepath: bool = True):
        """
        Method to create all required directories to run NLL.

        :param delete_basepath: Full pathname to save results from NonLinLoc
        """
        if os.path.isdir(self.nll_basepath) and delete_basepath is True:
            shutil.rmtree(self.nll_basepath)

        for dir_name in ["conf", "model", "obs", "time", "locs"]:
            if os.path.isdir(os.path.join(self.nll_basepath, dir_name)) is False:
                os.makedirs(os.path.join(self.nll_basepath, dir_name))

    def create_obs(self,
                   delete_file: bool = False):
        """
        Creates obs file for NLL.
        One file that contains all observations.

        :param delete_file: If True and .obs file exists, the existing .obs file will be deleted.
                            Default is False.
        """
        s = "{:<5} ? {:<3} ? {:<3} ? {:>25} GAU {:>6} -1 -1 -1 1\n"  # Format string

        # Check whether file exist
        if delete_file:
            try:
                open(os.path.join(self.nll_basepath, "obs", "events.obs"), "x")
            except FileExistsError:
                os.remove(os.path.join(self.nll_basepath, "obs", "events.obs"))  # Delete existing .obs file

        # Create one obs file that contains all events (including all available picks)
        with open(os.path.join(self.nll_basepath, "obs", "events.obs"), "w") as f_obs:
            for event in self.catalog.events:
                f_obs.write(f"# Origin time: {event.origins[0].time}\n")
                for pick in event.picks:
                    datetime_str = pick.time.datetime.strftime("%Y%m%d %H%M %S.%f")

                    # Read pick error uncertainty from picks in catalog. If not set, use default value from __init__.
                    uncertainty= pick.time_errors.uncertainty
                    if not uncertainty:
                        uncertainty = self.pick_uncertainty

                    #  Write each event to .obs file
                    f_obs.write(s.format(pick.waveform_id.station_code,
                                         "?",
                                         pick.phase_hint.upper(),
                                         datetime_str,
                                         uncertainty)
                                )
                f_obs.write("\n\n")  # Add to blank lines after each event in .obs file

    def create_config(self,
                      phase: str = "P",
                      signature: str = "LOCSIG AIS Phase Association",
                      comment: str = "LOCCOM AIS Phase Association"):
        """
        Creates config file for NonLinLoc.
        Perhaps further enhancements might be necessary, since this method only writes the required lines for the
        NonLinLoc config file.

        :param phase: Seismic phase (either P or S). Default is P, but both phase are required for NLL. So this method
                      runs twice
        :param signature: Signature string of NLL. Always starts with 'LOCSIG...'
        :param comment: Comment string of NLL. Always starts with 'LOCCOM...'
        """
        if phase.upper() not in ["P", "S"]:
            msg = "Parameter gtfile must be either 'S' or 'P' and not {}".format(phase)
            raise ValueError(msg)

        if signature.split()[0] != "LOCSIG" and len(signature.split()) <= 1:
            msg = "signature must start with 'LOCSIG' and a further signature."
            raise ValueError(msg)

        if comment.split()[0] != "LOCCOM" and len(comment.split()) <= 1:
            msg = "comment must start with 'LOCCOM' and a further comment."
            raise ValueError(msg)

        # Load velocity model
        self.__read_velocity_model()

        # Write config file
        with open(os.path.join(self.nll_basepath, "conf", f"{self.config_name}.conf"), "w") as f_conf:
            f_conf.write("# Generic control file statements\n")
            f_conf.write("CONTROL 2 51904\n")
            # Write earthquake location from phase associator
            #f_conf.write(f"TRANS SIMPLE {self.event.origins[0].latitude} {self.event.origins[0].longitude} 0.0\n")
            f_conf.write(f"TRANS SIMPLE {self.latorig} {self.longorig} 0.0\n")
            f_conf.write("# END of Generic control file statements\n")
            # Write Vel2Grid control file statements
            f_conf.write("# Vel2Grid control file statements\n")
            f_conf.write(f"VGOUT {os.path.join(self.nll_basepath, 'model', self.config_name)}\n")
            f_conf.write("VGTYPE P\n")
            f_conf.write("VGTYPE S\n")
            f_conf.write("VGGRID 2 2000 7000 0.0 0.0 0.0 0.5 0.5 0.5 SLOW_LEN\n")
            # f_conf.write("VGGRID 2 10000 7000 0.0 0.0 0.0 0.5 0.5 0.5 SLOW_LEN\n")  # XXX Probably this needs adjustment
            # Write velocity model
            for ind in range(len(self.df_vel_model)):
                f_conf.write(f"LAYER {self.df_vel_model['depth'][ind]} {self.df_vel_model['vp'][ind]} 0 "
                             f"{self.df_vel_model['vs'][ind]} 0 {self.df_vel_model['density'][ind]} 0\n")
            f_conf.write("# END of Vel2Grid control file statements\n")
            # Write Grid2Time control statements
            f_conf.write("# Grid2Time control file statements\n")
            f_conf.write(f"GTFILES {os.path.join(self.nll_basepath, 'model', self.config_name)} "
                         f"{os.path.join(self.nll_basepath, 'time', self.config_name)} {phase.upper()}\n")
            f_conf.write("GTMODE GRID2D ANGLES_YES\n\n")

            # Write all station_picks names and locations from json file to config
            for index in range(len(self.df_stations)):
                f_conf.write(f"GTSRCE {self.df_stations['id'][index].split('.')[1]} LATLON "
                             f"{self.df_stations['latitude'][index]} {self.df_stations['longitude'][index]} 0 0\n")

            f_conf.write("\nGT_PLFD 1.0e-3 0\n")
            f_conf.write("# END of Grid2Time control file statements\n")
            # Write NLL control file
            f_conf.write("# NLLoc control file statements\n")
            f_conf.write(f"{signature}\n")
            f_conf.write(f"{comment}\n")
            f_conf.write(f"LOCFILES {os.path.join(self.nll_basepath, 'obs', '*.obs')} NLLOC_OBS "
                         f"{os.path.join(self.nll_basepath, 'time', self.config_name)} "
                         f"{os.path.join(self.nll_basepath, 'locs', 'loc')}\n\n")
            f_conf.write("LOCPHASEID P\n")
            f_conf.write("LOCPHASEID S\n\n")
            f_conf.write("LOCGRID 6000 6000 1400 -200.0 -250.0 0.0 0.1 0.1 0.1 PROB_DENSITY SAVE\n\n")  # XXX
            f_conf.write("LOCSEARCH OCT 96 48 6 0.005 50000 10000 4 0\n\n")
            f_conf.write("LOCMETH EDT_OT_WT_ML 1.0e6 1 50 -1 -1.78 6 -1.0\n\n")
            f_conf.write("LOCQUAL2ERR 0.01 0.05 0.1 0.2 99999.9\n")
            f_conf.write("LOCANGLES ANGLES_NO 5\n")
            f_conf.write("LOCGAU 0.5 0.0\n")
            f_conf.write("LOCGAU2 0.01 0.05 2.0\n")
            f_conf.write("LOCHYPOUT SAVE_NLLOC_ALL\n\n")
            f_conf.write("# END of NLLoc control file statements\n")

    def __read_velocity_model(self):
        """
        Reading NLL velocity model as a pandas Dataframe object.
        The velocity model file should look like this:
            vp, vs, depth
            5, 3, 2
            5.5, 3.2, 5
            ...
        """
        if isinstance(self.vel_model, str):
            self.df_vel_model = pd.read_csv(self.vel_model)
        elif isinstance(self.vel_model, pd.DataFrame):
            self.df_vel_model = self.vel_model

    def run_vel2grid(self):
        """
        Running NLL Vel2Grid.
        """
        execute = os.path.join(self.nll_executables, 'Vel2Grid')
        os.system(f"{execute} {self.conf_file}")

    def run_grid2time(self):
        """
        Running NLL Grid2Time for phases P and S.
        """
        execute = os.path.join(self.nll_executables, 'Grid2Time')
        for phase in ["P", "S"]:
            if phase == "S":
                self.create_config(phase="S")
            os.system(f"{execute} {self.conf_file}")

    def localise(self):
        """
        Running NLL NLLoc for localisation of events.
        """
        execute = os.path.join(self.nll_executables, 'NLLoc')
        os.system(f"{execute} {self.conf_file}")

    def run_nll(self,
                delete_dirs: Optional[list] = None):
        """
        Run all NLL method for a complete localisation. i.e.:
        1. Vel2Grid
        2. Grid2Time for P and S
        3. NLLoc

        Since NLL model files can require much memory, these files can be deleted after one run.
        However, if these files are not deleted, then localisation of events is much faster, since these files are
        only created once.

        :param delete_dirs: List of directories that will be deleted. Per default this value is None.
        """
        self.run_vel2grid()
        self.run_grid2time()
        self.localise()

        # Delete travel time files in directory time
        if delete_dirs:
            for dir_name in delete_dirs:
                self.delete_dir_content(directory=dir_name)

    def delete_dir_content(self,
                           directory: str):
        """
        Method to delete the whole content of a directory using Python's os module.

        :param directory: Directory name to delete.
        """
        files = glob.glob(os.path.join(self.nll_basepath, directory, "*"))
        for filename in files:
            os.remove(filename)


def nll_object(nll):
    """
    Deletes 'time' and 'model' directories to save memory storage. However, if these files are deleted, they have to be
    recreated for each run of NonLinLoc which is time-intensive for large seismic networks.
    """
    nll.run_nll(delete_dirs=["time", "model"])


def update_events_from_nll(station_df: Union[str, pd.DataFrame],
                           nll_basepath: str,
                           depth_filter: float = 10000) -> list:
    """
    Reads created event files from NonLinLoc and returns a list that contains all found events.

    :param station_df: csv-filename or pandas dataframe that contains station information
    :param nll_basepath: Full pathname to save results from NonLinLoc
    :param depth_filter: Removes event with depth uncertainty greater than the given threshold. Default are 10,000 m.
    """
    # Read all .hyp files from NonLinLoc
    hyp_files = glob.glob(os.path.join(nll_basepath, 'locs', 'loc.*.grid0.loc.hyp'))

    # Remove sum file from .hyp files and read each fname as a single event and append to event_lst
    event_lst = []
    for fname in hyp_files:
        if os.path.split(fname)[1].split(".")[1] != "sum":
            nll_event = read_events(pathname_or_url=fname,
                                    format="NLLOC_HYP")[0]
            # Apply depth filter, i.e. if depth uncertainty of nll_event is greater than a threshold, the event
            # will not be used further. # Default value (5000 m) is taken from S-EqT Paper
            if depth_filter:
                depth_uncertainty = nll_event.origins[0].depth_errors.uncertainty
                if depth_uncertainty <= depth_filter:
                    event_lst.append(nll_event)
            else:
                event_lst.append(nll_event)

    # Since NonLinLoc only works with station_picks names and not with station_picks and network names, the network needs to
    # be added to each single pick for each event in event_lst
    # Load station_json as pandas dataframe
    if isinstance(station_df, pd.DataFrame):
        station_df = station_df
    else:
        station_df = load_stations(station_json=station_df)

    # Loop over each event and pick and add network and location to waveform_id
    for event in event_lst:
        for pick in event.picks:
            for id in station_df["id"]:
                if id.split(".")[1] == pick.waveform_id.station_code:
                    network, station, location = id.split(".")
                    break
            pick.waveform_id.network_code = network
            pick.waveform_id.location_code = location

    return event_lst


def check_nll_time(station_json: str,
                   nll_basepath: str) -> bool:
    """
    Check whether for each station_picks all required files in time directory exist.
    If all files exist, function returns True, otherwise False
    :param station_json:
    :param nll_basepath:
    """
    if isinstance(station_json, pd.DataFrame):
        df_stations = station_json
    else:
        df_stations = load_stations(station_json=station_json)

    for station_id in df_stations["id"]:
        name = station_id.split(".")[1]
        for phase in ["P", "S"]:
            for ftype in ["angle.buf", "angle.hdr", "time.buf", "time.hdr"]:
                try:
                    filename = glob.glob(os.path.join(nll_basepath, "time", f"*{phase}.{name}.{ftype}"))[0]
                except IndexError:
                    return False
                if os.path.isfile(filename) is False:
                    return False

    return True


def nll_wrapper(catalog: obspy.Catalog,
                station_df: Union[str, pd.DataFrame],
                nll_basepath: str = "../NonLinLoc",
                vel_model: str = "../metadata/velocity.nll",
                nll_executable: str = "/work/software/nlloc/7.00.16/src/bin") -> obspy.Catalog:
    """
    Wrapper to relocalise an obspy Catalog with NonLinLoc. The obspy catalog can be created, for example, by a phase
    association algorithm. Picks must be included in the catalog. Additionally, stations are given either as pandas
    Dataframe or the Dataframe will be created from the give filename. The station_df must contain the following rows,
    including the header: trace_id (i.e. network_name.station_name.location), latitude, longitude, elevation_m).
    The velocity model file should look like this (vp, vs in km/s and depth in km):
        vp, vs, depth
        5, 3, 2
        5.5, 3.2, 5
        ......

    Note, if model and time files exist in, NonLinLoc will use these existing files for relocalisation. Otherwise,
    these files will be created.

    :param catalog: Catalog to relocalise data
    :param station_df: csv-filename or pandas dataframe that contains station information
    :param nll_basepath: Full pathname to save results from NonLinLoc
    :param vel_model: Full pathname of velocity model
    :param nll_executable: Executable of NonLinLoc

    :return: obspy Catalog with relocalised seismic events.
    """
    # Check files in NLL basepath directory
    nll_time_files_exist = check_nll_time(station_json=station_df,
                                          nll_basepath=nll_basepath)
    if nll_time_files_exist is True:
        create_nll_files = False
    else:
        create_nll_files = True

    # Create NonLinLoc instance from the NLL Class
    nll = Event2NLL(catalog=catalog,
                    nll_basepath=nll_basepath,
                    vel_model=vel_model,
                    stations_df=station_df,
                    nll_executable=nll_executable,
                    create_files=create_nll_files)

    if nll_time_files_exist is True:
        nll.create_obs(delete_file=True)  # Create NLL .obs file
        nll.delete_dir_content(directory=os.path.join(nll_basepath, "locs"))  # Delete locs directory
        nll.localise()  # Localise events
    else:
        nll.run_nll()  # Do full run of NLL, i.e. Vel2Grid, Grid2Time and NLLoc

    # Update each localisation of each event
    # Catalog is also filtered with depth uncertainty
    all_events = update_events_from_nll(station_df=station_df,
                                        nll_basepath=nll_basepath)

    # Export full catalog
    full_catalog = obspy.Catalog(events=sort_events(all_events))

    return full_catalog
