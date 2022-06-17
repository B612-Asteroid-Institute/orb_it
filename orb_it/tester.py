import os
import logging
import numpy as np
import pandas as pd
import astropy
from astropy.time import Time
from astropy import units as u
from astroquery.jplhorizons import Horizons
from shutil import copy
from .raiden import Orbits

logging.basicConfig(handlers=[logging.FileHandler('testing.log'),logging.StreamHandler()],level=logging.INFO,format='%(asctime)s, %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

def error_state(orbit, observatory_code, backend, astrometric_error, full_output,ephemeris=None):
    result = pd.DataFrame()
#     result["num_obs_fit"] = residuals["incl"].sum().astype(int)
    result['epoch [mjd]'] = [np.nan]
    result.insert(0,'orbit_id',orbit.ids[0])
    result.insert(1,'integrator',backend.name)
    result["delta epoch [mjd]"] = [np.nan]
    result["delta r [km]"] = [np.nan]
    result["delta v [m/s]"] = [np.nan]
    result["delta x [km]"] = [np.nan]
    result["delta y [km]"] = [np.nan]
    result["delta z [km]"] = [np.nan]
    result["delta vx [m/s]"] = [np.nan]
    result["delta vy [m/s]"] = [np.nan]
    result["delta vz [m/s]"] = [np.nan]
    if full_output:
        result["rms delta ra [arcsec]"] = [np.nan]
        result["rms delta dec [arcsec]"] = [np.nan]
        result["rms delta time [seconds]"] =[np.nan]
    
        result["covariance"] = [np.nan]
    
    
    result.insert(2, "observatory_code",  observatory_code)
    result.insert(3, "arc_length [days]", np.nan)
    result.insert(6, "astrometric_error [mas]", astrometric_error)
    if ephemeris is not None:
        result.insert(4, "num_obs", len(ephemeris))
    else:
        result.insert(4, "num_obs", np.nan)
    return result

def doTest(orbit, observatory_code, dts, backend, astrometric_error=None,  full_output=False,out=None):
    t0 = orbit.epochs[0]
    if isinstance(dts, astropy.time.core.Time):
        #assert t0 == dts[0], "The first time in the dts should be the observed time of your initial orbits"
        observation_times = dts
    elif isinstance(dts, (np.ndarray,list)):
        observation_times = t0 + dts
    if astrometric_error is None:
        astrometric_error = 0
    if out is None:
        outd=None
    else:
        whq1='_'.join(orbit.ids[0].split(' '))
        whq2 = whq1.split('/')
        nam='_'.join(whq2)
        tm=Time.now().utc.strftime('%Y%m%d_%H%M')
        outd=os.path.join(out,f"{nam}",backend.name,f"{dts.max():0.0f}_{astrometric_error:0.0f}_{tm}")
        #outd=os.path.join(out,f"{nam}/{dts.max():0.0f}_{astrometric_error:0.0f}mas_{Time.now().utc.value.isoformat(timespec='seconds')}")
    MAS_TO_DEG = 2.777777777777778e-07
    DEG_TO_MAS = 1/MAS_TO_DEG
    obs = {observatory_code:observation_times}
    try:
        if full_output:
            ephemeris,ret1= backend._generateEphemeris(orbits=orbit,observers =obs, out_dir=outd,full_output=full_output)
        else:
            ephemeris= backend._generateEphemeris(orbits=orbit,observers =obs, out_dir=outd)
    except:
        logging.error(f"{orbit.ids[0]} failed to run _generateEphemeris using {backend.name}")
        #print(f"Error in {backend.name} _generateEphemeris")
        return error_state(orbit, observatory_code, backend, astrometric_error, full_output)

    ephemeris["RA_sigma_deg"] = astrometric_error*MAS_TO_DEG
    ephemeris["Dec_sigma_deg"] = astrometric_error*MAS_TO_DEG
    ephemeris["RA_sigma_mas"] = astrometric_error
    ephemeris["Dec_sigma_mas"] = astrometric_error
    ephemeris["mjd_sigma_seconds"] = [0.01 for i in range(len(ephemeris))]
    ephemeris["obs_id"] = [f"o{i:07d}" for i in range(len(ephemeris))]
    #ephemeris.rename(columns={"mjd_utc" : "mjd", "orbit_id" : "obj_id"}, inplace=True)

    # Add a simple astrometric errors to both RA and Dec
    ephemeris.loc[:, "RA_deg"] += np.random.normal(loc=0, scale=astrometric_error*MAS_TO_DEG, size=len(ephemeris)) / np.cos(np.radians(ephemeris['Dec_deg'].values))
    ephemeris.loc[:, "Dec_deg"] += np.random.normal(loc=0, scale=astrometric_error*MAS_TO_DEG, size=len(ephemeris))
    try:
        if full_output:
            od_orbit_df, residuals,ret2 = backend._orbitDetermination(ephemeris, out_dir=outd,full_output=full_output)
        else:
            od_orbit_df = backend._orbitDetermination(ephemeris, out_dir=outd)
    
        od_orbit = Orbits(od_orbit_df)
    except:
        logging.error(f"{orbit.ids[0]} failed to run _orbitDetermination using {backend.name}")
        #print(f"Error in {backend.name} _orbitDetermination")
        return error_state(orbit, observatory_code, backend, astrometric_error, full_output, ephemeris)

    try:
        if full_output:
            prop_orbit,ret3 = backend._propagateOrbits(orbit, observation_times[-1:], out_dir=outd,full_output=full_output)
        else:
            prop_orbit = backend._propagateOrbits(orbit, od_orbit.epochs, out_dir=outd)
    except:
        logging.error(f"{orbit.ids[0]} failed to run _propagateOrbits using {backend.name}")
        #print(f"Error in {backend.name} _propagateOrbits")
        return error_state(orbit, observatory_code, backend, astrometric_error, full_output, ephemeris)

    delta_state = od_orbit_df[["x", "y", "z", "vx", "vy", "vz"]].values - prop_orbit[["x", "y", "z", "vx", "vy", "vz"]].values
    
    delta_r = np.linalg.norm(delta_state[:,0:3])
    delta_v = np.linalg.norm(delta_state[:,3:6])
    
    # os.rmtree, rmdir add parameter check and dir deletion
    
    result = pd.DataFrame()
    
    result['orbit_id'] = od_orbit_df.orbit_id.values
    result['integrator'] = backend.name
#     result["num_obs_fit"] = residuals["incl"].sum().astype(int)
    result['epoch [mjd]'] = prop_orbit.mjd_tdb.values
    result["delta epoch [mjd]"] = od_orbit.epochs.tdb.mjd - prop_orbit["mjd_tdb"].values
    result["delta r [km]"] = (delta_r * u.AU).to(u.km).value
    result["delta v [m/s]"] = (delta_v * u.AU / u.d).to(u.m / u.s).value 
    result["delta x [km]"] = (delta_state[:,0] * u.AU).to(u.km).value
    result["delta y [km]"] = (delta_state[:,1] * u.AU).to(u.km).value
    result["delta z [km]"] = (delta_state[:,2] * u.AU).to(u.km).value
    result["delta vx [m/s]"] = (delta_state[:,3] * u.AU / u.d).to(u.m / u.s).value
    result["delta vy [m/s]"] = (delta_state[:,4] * u.AU / u.d).to(u.m / u.s).value
    result["delta vz [m/s]"] = (delta_state[:,5] * u.AU / u.d).to(u.m / u.s).value
    if full_output:
        result["rms delta ra [arcsec]"] = np.sqrt(np.mean(residuals["dRA"].values**2))
        result["rms delta dec [arcsec]"] = np.sqrt(np.mean(residuals["dDec"].values**2))
        result["rms delta time [seconds]"] = np.sqrt(np.mean(residuals["dTime"].values**2))
    
        result["covariance"] = od_orbit_df["covariance"].values
    
    arc_length = observation_times.utc.mjd.max() - observation_times.utc.mjd.min()
    
    result.insert(2, "observatory_code",  observatory_code)
    result.insert(3, "arc_length [days]", arc_length)
    result.insert(6, "astrometric_error [mas]", astrometric_error)
    result.insert(4, "num_obs", len(ephemeris))
    # if out is not None and full_output:
    #     bashgen(outd,astrometric_error,observation_times[-1:].utc.jd,observatory_code)
    
    return result

def runTest(orbits, observatory_code, dts, backend, astrometric_error=None, full_output=False, out=None):
    '''
    Runs the end-to-end test for a list of orbits. Generates ephemeris, 
    conducts orbit determination, and propagates the initial orbits to the final 
    time and returns the comparison.

    Parameters
    ----------
    orbits : orbit object `~orb_it.raiden.Orbits`
        Orbits to be tested.
    observatory_code : str or list of strings
        MPC Observatory code(s) for the observatory to be tested. (500 for Geocenter)
    dts : array of floats, 2D array of floats, or list of `~astropy.time.core.Time` objects with the same length as `orbits`
        List of observation times after the initial time to test the orbits over. Measured in days.
        NOTE: Anything passed to this parameter must have values in ASCENDING ORDER.
    backend : backend object `~orb_it.backend.backend.BACKEND`
        Integrator backend to be used for the test.
    astrometric_error : float or list of floats, optional
        Astrometric error to be added to the generated observations. 
        If None, no astrometric error is added. Units are milliarcseconds.
    full_output : bool, optional
        If True, returns the full output of the test including residuals and covariance matrices.
        DOES NOT WORK WITH PYOORB CURRENTLY.
    out : str, optional
        Path to the output directory for saving necessary files for this test to be
        run independently. This includes configuration files, generated files, and bash scripts. 
        It has the file structure '{out}/{orbit_id}/{days propagated}days_{error}mas_{timestamp}'. 
        If None, no files are saved.
        DOES NOT WORK WITH PYOORB CURRENTLY.
    
    Returns
    -------
    result : `~pandas.DataFrame`
        DataFrame with the results of the test.
        Has the following columns:
        orbit_id : str
            Orbit ID of the object being tested.
        integrator : str
            Name of the integrator backend used for test.
        observatory_code : str
            MPC Observatory code for the observatory to be tested. (500 for Geocenter)
        arc_length [days] : float
            Length of time in days over which the orbits are tested.
        num_obs : int
            Number of observations generated in the test.
        epoch [mjd] : float
            Time at the end of the arc in mjd (t0 + final dt day).
        astrometric_error [mas] : float
            Astrometric error added to the observations.
        delta epoch [mjd] : float
            Difference between the final epochs of the orbit determination and the propagated orbit.
        delta r [km] : float
            The absolute distance between the orbit determination result and the propagated orbit.
            Measured in kilometers.
        delta v [m/s] : float
            The absolute velocity difference between the orbit determination result and the propagated orbit.
            Measured in meters per second.
        delta x [km] : float
            The x position difference between the orbit determination result and the propagated orbit.
            Measured in kilometers.
        delta y [km] : float
            The y position difference between the orbit determination result and the propagated orbit.
            Measured in kilometers.
        delta z [km] : float
            The z position difference between the orbit determination result and the propagated orbit.
            Measured in kilometers.
        delta vx [m/s] : float
            The x velocity difference between the orbit determination result and the propagated orbit.
            Measured in meters per second.
        delta vy [m/s] : float
            The y velocity difference between the orbit determination result and the propagated orbit.
            Measured in meters per second.
        delta vz [m/s] : float
            The z velocity difference between the orbit determination result and the propagated orbit.
            Measured in meters per second.
        rms delta ra [arcsec] : float, choose `full_output=True` to get this value
            The root mean square difference of the Right Acension of the Residuals from Orbit Determination.
            Measured in arcseconds.
        rms delta dec [arcsec] : float, choose `full_output=True` to get this value
            The root mean square difference of the Declination of the Residuals from Orbit Determination.
            Measured in arcseconds.
        rms delta time [seconds] : float, choose `full_output=True` to get this value
            The root mean square difference of the Time of the Residuals from Orbit Determination.
            Measured in seconds.
        covariance : 2D array converted to 3D array, choose `full_output=True` to get this value
            Covariance matrix of the orbit determination result.
    '''
    if isinstance(observatory_code,(list,np.ndarray)) and isinstance(observatory_code[0],str):
        pass
    elif isinstance(observatory_code,str):
        observatory_code = [observatory_code]
    
    try:
        astrometric_error[0]
        errors = astrometric_error
    except:
        errors = [astrometric_error]
    
    #1-D cases
    #array of ints np.ndarray: [1,2,3,...]
    if isinstance(dts, (list,np.ndarray)) and isinstance(dts[0], (int,float,np.int64,np.float64)):
        dts = [dts]
    #list in time object: <Time object: [t1,t2,...]>
    elif isinstance(dts, astropy.time.core.Time) and isinstance(dts[0], astropy.time.core.Time):
        dts = [dts]
    #2-D cases
    #2d array of ints: [[1,2,3,...],[...]]
    elif isinstance(dts, (list,np.ndarray)) and isinstance(dts[0], (list,np.ndarray)):
        pass
    #2d array of time objects: list([<Time object: [t1,t2,...]>,<Time object: [t1,t2,...]>,...])
    elif isinstance(dts, (list,np.ndarray)) and isinstance(dts[0], astropy.time.core.Time) and len(dts) == orbits.num_orbits:
        
        results_i = []
        for i in range(orbits.num_orbits):
        #this will iterate over unique dts for each orbit
            for k in range(len(errors)):
                for l in range(len(observatory_code)):
                    result = doTest(orbits[i],observatory_code[l],dts[i],astrometric_error=errors[k],out=out, full_output=full_output, backend=backend)
                    results_i.append(result)
        results = pd.concat(
            results_i,
            ignore_index=True
        )
        return results
    else:
        raise ValueError('dts must be a list of floats, array of time objects, 2D array of floats, or 2D array of time objects corresponding to each orbit being tested.')
    #try:
        # check for 2d array
      #  dts[0][0]
    #except:
    #    dts = [dts]
    results_i = []
    for i in range(orbits.num_orbits):
        for j in range(len(dts)):
            for k in range(len(errors)):
                for l in range(len(observatory_code)):
                    result = doTest(orbits[i],observatory_code[l],dts[j],astrometric_error=errors[k],out=out, full_output=full_output, backend=backend)
                    results_i.append(result)
    results = pd.concat(
        results_i,
        ignore_index=True
    )
    return results