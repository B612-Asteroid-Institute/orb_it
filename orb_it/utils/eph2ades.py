import numpy as np
import pandas as pd
from astropy.time import Time
import os

__all__ = [
    "writeADESHeader",
    "writeToADES"
]
class ADES:
    def __init__(self,eph,ades,error,t_error,obs_code):
        self.eph = eph
        self.ades = ades
        self.error = error
        self.t_error = t_error
        self.observatory_code = obs_code

    def ephMod(self):
        if (self.observatory_code == "500"):
                columns = [
                        "jd_utc", "RA_deg", "Dec_deg", "delta_au", "r_au",
                        "elong", "ph_ang", "mag",
                        "'/hr", "PA", "rvel",
                        "lon", "lat", "altitude_km"
                    ]
        else:
                    columns = [
                        "jd_utc", "RA_deg", "Dec_deg", "delta_au", "r_au",
                        "elong", "ph_ang", "mag",
                        "RA_'/hr", "dec_'/hr",
                        "alt",  "az",  "rvel",
                        "lon", "lat", "altitude_km"
                    ]
        MAS_TO_DEG = 2.777777777777778e-07
        DEG_TO_MAS = 1/MAS_TO_DEG
        ephemeris=pd.read_csv(
                            self.eph,
                            header=0,
                            delim_whitespace=True,
                            names=columns,
                            float_precision="round_trip"

                        )
        times = Time(
            ephemeris["jd_utc"].values,
            format="jd",
            scale="utc"
        )
        ephemeris["mjd_utc"] = times.utc.mjd
        ephemeris.drop(
            columns=["jd_utc"],
            inplace=True
        )
        ephemeris['observatory_code'] = np.full(len(ephemeris["RA_deg"]),self.observatory_code)
        astrometric_error=self.error
        ephemeris["RA_sigma_deg"] = astrometric_error*MAS_TO_DEG
        ephemeris["Dec_sigma_deg"] = astrometric_error*MAS_TO_DEG
        ephemeris["mjd_sigma_seconds"] = [self.t_error for i in range(len(ephemeris))]
        ephemeris["obs_id"] = [f"o{i:07d}" for i in range(len(ephemeris))]
        #ephemeris.rename(columns={"mjd_utc" : "mjd", "orbit_id" : "obj_id"}, inplace=True)

        # Add a simple astrometric errors to both RA and Dec
        ephemeris.loc[:, "RA_deg"] += np.random.normal(loc=0, scale=astrometric_error*MAS_TO_DEG, size=len(ephemeris)) / np.cos(np.radians(ephemeris['Dec_deg'].values))
        ephemeris.loc[:, "Dec_deg"] += np.random.normal(loc=0, scale=astrometric_error*MAS_TO_DEG, size=len(ephemeris))
        self.observationsToADES(ephemeris)
        return

    def observationsToADES(self, observations):
    
        #if not os.path.exists(out_dir):
        #    os.mkdir(out_dir)
        
        _observations = observations.copy()
        _observations.rename(
            columns={
                "mjd_utc" : "mjd",
                "RA_deg" : "ra",
                "Dec_deg" : "dec",
                "RA_sigma_deg" : "rmsRA",
                "Dec_sigma_deg" : "rmsDec",
                "mag_sigma" : "rmsMag",
                "mjd_sigma_seconds" : "rmsTime",
                "filter" : "band",
                "observatory_code" : "stn",
            },
            inplace=True
        )
        
        # ADES expects rms in arcseconds not degrees
        _observations["rmsRA"] = _observations["rmsRA"].values * 3600
        _observations["rmsDec"] = _observations["rmsDec"].values * 3600
        _observations.sort_values(
            by=["mjd", "stn"],
            inplace=True
        )

        #id_present = False
        #id_col = None
        #if "permID" in _observations.columns.values:
        #    id_present = True
        #    id_col = "permID"
        #if "provID" in _observations.columns.values:
        #    id_present = True
        #    id_col = "provID"
        #if "trkSub" in _observations.columns.values:
        #    id_present = True
        #    id_col = "trkSub"
        #if not id_present:
        #    if "orbit_id" in _observations.columns:
        #        _observations["trkSub"] = _observations["orbit_id"]
        #    else:
        #        _observations["trkSub"] = _observations["obj_id"]
        #    id_col = "trkSub"
        id_col="trkSub"
        _observations["trkSub"] = np.full(len(_observations['ra']),'o00000000')
        if "mag" not in _observations.columns:
            _observations.loc[:, "mag"] = 20.0
        if "band" not in _observations.columns:
            _observations.loc[:, "band"] = "V"
        if "mode" not in _observations.columns:
            _observations.loc[:, "mode"] = "CCD"
        if "astCat" not in _observations.columns:
            _observations.loc[:, "astCat"] = "None"


        for i, orbit_id in enumerate(_observations[id_col].unique()):
            # If you give fo a string for a numbered object it will typically append brackets
            # automatically which makes retrieving the object's orbit a little more tedious so by making sure
            # the object ID is not numeric we can work around that.
            orbit_id_i = f"o{i:08d}"

            observations_file = self.ades
            #out_dir = os.path.join(out_dir, "od_{}".format(orbit_id_i))

            mask = _observations[id_col].isin([orbit_id])
            object_observations = _observations[mask].copy()
            object_observations.loc[:, id_col] = orbit_id_i
            object_observations.reset_index(inplace=True, drop=True)

            self.writeToADES(
                object_observations,
                observations_file,
                mjd_scale="utc",
                seconds_precision=9,
                columns_precision={
                "ra" : 16,
                "dec" : 16,
                "mag" : 2,
                "rmsMag" : 2}
            )
                
        return

    

    def writeADESHeader(
        self,
            observatory_code,
            submitter,
            telescope_design,
            telescope_aperture,
            telescope_detector,
            observers,
            measurers,
            observatory_name=None,
            submitter_institution=None,
            telescope_name=None,
            telescope_fratio=None,
            comment=None
        ):
        """
        Write the ADES PSV headers.

        Parameters
        ----------
        observatory_code : str
            MPC-assigned observatory code
        submitter : str
            Submitter's name.
        telescope_design : str
            Telescope's design, eg. Reflector.
        telescope_aperture : str
            Telescope's primary aperture in meters.
        telescope_detector : str
            Telescope's detector, eg. CCD.
        observers : list of str
            First initial and last name (J. Smith) of each of the observers.
        measurers : list of str
            First initial and last name (J. Smith) of each of the measurers.
        observatory_name : str, optional
            Observatory's name.
        submitter_insitution : str, optional
            Name of submitter's institution.
        telescope_name : str, optional
            Telescope's name.
        telescope_fratio : str, optional
            Telescope's focal ratio.
        comment : str
            Additional comment to add to the ADES header.


        Returns
        -------
        list : str
            A list of each line in the ADES header.
        """
        # Start header with version number
        header = [
            "# version=2017",
        ]

        # Add observatory [required]
        header += ["# observatory"]
        header += [f"! mpcCode {observatory_code}"]
        if observatory_name is not None:
            header += [f"! name {observatory_name}"]

        # Add submitter [required]
        header += ["# submitter"]
        header += [f"! name {submitter}"]

        if submitter_institution is not None:
            header += ["! institution {}".format(submitter_institution)]

        # Add telescope details [required]
        header += ["# telescope"]
        if telescope_name is not None:
            header += [f"! name {telescope_name}"]
        header += [f"! design {telescope_design}"]
        header += [f"! aperture {telescope_aperture}"]
        header += [f"! detector {telescope_detector}"]
        if telescope_fratio is not None:
            header += [f"! fRatio {telescope_fratio}"]

        # Add observer details
        header += ["# observers"]
        if type(observers) is not list:
            err = (
                "observers should be a list of strings."
            )
            raise ValueError(err)
        for name in observers:
            header += [f"! name {name}"]

        # Add measurer details
        header += ["# measurers"]
        if type(measurers) is not list:
            err = (
                "measurers should be a list of strings."
            )
            raise ValueError(err)
        for name in measurers:
            header += [f"! name {name}"]

        # Add comment
        if comment is not None:
            header += ["# comment"]
            header += ["! line {}".format(comment)]

        header = [i + "\n" for i in header]
        return header

    def writeToADES(
        self,
            observations,
            file_out,
            mjd_scale="utc",
            seconds_precision=9,
            columns_precision={
                "ra" : 16,
                "dec" : 16,
                "mag" : 2,
                "rmsMag" : 2,
            },
            observatory_code="I11",
            submitter="D. iRAC",
            telescope_design="Reflector",
            telescope_aperture="8.4",
            telescope_detector= "CCD",
            observers=["D. iRAC"],
            measurers=["D. iRAC"],
            observatory_name="Vera C. Rubin Observatory",
            submitter_institution=None,
            telescope_name=None,
            telescope_fratio=None,
            comment=None
        ):
        """
        Save observations to a MPC-submittable ADES psv file.

        Parameters
        ----------
        observations : `~pandas.DataFrame`
            Dataframe containing observations.
        file_out : str
            Path and name to save out
        mjd_scale : str, optional
            Time scale of MJD observation times
        seconds_precision : int, optional
            Number of decimal places of precision on the measurement
            of seconds for the observation times. The ADES format can handle higher
            than ms precision if the observations warrant such accuracy. 0.1 ms precision
            would be expressed as 4 while ms precision would be expressed as 3.
        columns_precision : dict, optional
            Dictionary with column names as keys and the precision (in decimals) to which
            they should be printed in the ADES file.
        observatory_code : str, optional
            MPC-assigned observatory code
        submitter : str, optional
            Submitter's name.
        telescope_design : str, optional
            Telescope's design, eg. Reflector.
        telescope_aperture : str, optional
            Telescope's primary aperture in meters.
        telescope_detector : str, optional
            Telescope's detector, eg. CCD.
        observers : list of str, optional
            First initial and last name (J. Smith) of each of the observers.
        measurers : list of str, optional
            First initial and last name (J. Smith) of each of the measurers.
        telescope_name : str, optional
            Telescope's name.
        telescope_fratio : str, optional
            Telescope's focal ratio.
        comment : str, optional
            Additional comment to add to the ADES header.

        Returns
        -------
        list : str
            A list of each line in the ADES header.
        """
        header = self.writeADESHeader(
            observatory_code,
            submitter,
            telescope_design,
            telescope_aperture,
            telescope_detector,
            observers,
            measurers,
            observatory_name=observatory_name,
            submitter_institution=submitter_institution,
            telescope_name=telescope_name,
            telescope_fratio=telescope_fratio,
            comment=comment
        )

        # Format columns from observations into PSV format
        ades = {}

        id_present = False
        if "permID" in observations.columns.values:
            ades["permID"] = observations["permID"].values
            id_present = True
        if "provID" in observations.columns.values:
            ades["provID"] = observations["provID"].values
            id_present = True
        if "trkSub" in observations.columns.values:
            ades["trkSub"] = observations["trkSub"].values
            id_present = True

        if not id_present:
            err = (
                "At least one of permID, provID, or trkSub should\n"
                "be present in observations."
            )
            raise ValueError(err)

        observation_times = Time(
            observations["mjd"].values,
            format="mjd",
            scale=mjd_scale,
            precision=seconds_precision
        )
        ades["obsTime"] = np.array([i + "Z" for i in observation_times.utc.isot])
        ades["ra"] = observations["ra"].values
        ades["dec"] = observations["dec"].values

        if "rmsRA" in observations.columns.values:
            ades["rmsRA"] = observations["rmsRA"].values
        if "rmsDec" in observations.columns.values:
            ades["rmsDec"] = observations["rmsDec"].values

        ades["mag"] = observations["mag"].values
        if "rmsMag" in observations.columns.values:
            ades["rmsMag"] = observations["rmsMag"].values
        if "rmsTime" in observations.columns.values:
            ades["rmsTime"] = observations["rmsTime"].values
        if "uncTime" in observations.columns.values:
            ades["uncTime"] = observations["uncTime"].values
        ades["band"] = observations["band"].values
        ades["stn"] = observations["stn"].values
        ades["mode"] = observations["mode"].values
        ades["astCat"] = observations["astCat"].values

        if "remarks" in observations.columns.values:
            ades["remarks"] = observations["remarks"].values

        for col in columns_precision:
            if col in ades.keys():
                prec_col = columns_precision[col]
                ades[col] = [f"{i:.{prec_col}f}" for i in ades[col]]

        # Create dataframe with formated data entries
        ades = pd.DataFrame(ades)
        col_header = "|".join(ades.columns) + "\n"

        with open(file_out, "w") as f:
            f.write("".join(header))
            f.write(col_header)

        ades = ades.replace(
            np.nan,
            " ",
            regex=True
        )

        #reduced_precision_cols = ["rmsMag", "uncTime", "rmsTime"]
        reduced_precision_cols = []
        for col in reduced_precision_cols:
            if col in ades.columns:
                ades[col] = ades[col].map(lambda x: '{0:.3f}'.format(x))

        ades.to_csv(
            file_out,
            sep="|",
            header=False,
            index=False,
            mode="a",
            float_format='%.16f'
        )
        return

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('ephemeris_path', type=str, help='Path to ephemeris file')
    parser.add_argument('ades_path', type=str, help='Path to new ADES file')
    parser.add_argument('--astrometric_error', type=float, default=0.0, help='Astrometric error')
    parser.add_argument('--time_error', type=float, default=0.01, help='Time error')
    parser.add_argument('--observatory_code', type=str, default='500', help='Observatory code')
    args = parser.parse_args()
    e2a = ADES(args.ephemeris_path, args.ades_path, args.astrometric_error, args.time_error, args.observatory_code)
    e2a.ephMod()

