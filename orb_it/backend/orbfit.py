from astropy.time import Time
import numpy as np
import pandas as pd
import tempfile as tf
import os
from astropy.coordinates import Angle
import astropy.units as u
import subprocess

from .backend import Backend
from ..utils.ic_utils import *

# As of right now, orbfit_path should only lead to HOME/orbfit.
ORBFIT_CONFIG = {'orbfit_path': os.path.join(os.environ["HOME"],'orbfit/'),
                 'arc_limit': 31
}

class ORBFIT(Backend):
    def __init__(self,  **kwargs):
        allowed_kwargs = ORBFIT_CONFIG.keys()
        for k in kwargs:
            if k not in allowed_kwargs:
                raise ValueError()

        # If an allowed kwarg is missing, add the
        # default
        for k in allowed_kwargs:
            if k not in kwargs:
                kwargs[k] = ORBFIT_CONFIG[k]


        super().__init__(name="OrbFit", **kwargs)

        self.setup()
        return

    def setup(self):
        if not os.path.exists(self.orbfit_path):
            FileNotFoundError(f'orbfit not found in {self.orbfit_path}')
        # Add extra copying things for the AST files and the jpleph file
        
        return

    def _generateEphemeris(self, orbits, observers, out_dir=None):
        '''
        Generate ephemerides for each orbit and observer.

        Parameters
        ----------
        orbits : `~orb_it.raiden.Orbits`
            Orbits to propagate.
        observers : dict
            A dictionary with observatory codes as keys and observation_times (`~astropy.time.core.Time`) as values.
        out_dir : str, optional
            Save input and output files to this directory. Will create
            a sub directory called ephemeris inside this directory.

        Returns
        -------
        ephemeris : `~pandas.DataFrame`
            Ephemerides for each orbit and observer.
        '''
        ephemeris_dfs = []
        for observatory_code, observation_times in observers.items():
            for i in range(orbits.num_orbits):
                orbit_id_path = "_".join(orbits.ids[i].astype(str).split(" "))
                orbit_id_path = "_".join(orbit_id_path.split("/"))
                orbit_id_i = orbits.ids[i]
                start = observation_times[0]
                stop = observation_times[-1]
                step = observation_times[1].mjd - observation_times[0].mjd
                with tf.TemporaryDirectory() as tdir:
                    eq1file(orbits[i], tdir=tdir)
                    fopfile(tdir=tdir,obsc=observatory_code)
                    ephHelp(start, stop, step, observatory_code, tdir=tdir)
                    home = os.environ["HOME"]
                    subprocess.call(f'{home}/orbfit/src/fitobs/fitobs.x < ast.inp',cwd=tdir,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL) 
                    f1 = open(f'{tdir}/temp.eph').read().split('\n')
                    for i in range(len(f1)):
                        if '=' in f1[i]:
                            d1 = f1[i+1:]
                            break
                    a1=[]
                    for i in range(len(d1)):
                        if "." in d1[i]:
                            a1.append(d1[i].split())
                    df = pd.DataFrame(a1)
                    
                df['orbit_id'] = orbit_id_i
                df['mjd_utc'] = df[4].astype(float)
                df["RA_deg"] = Angle((df[5].astype(float), df[6].astype(float), df[7].astype(float)), unit = 'hourangle').degree
                df["Dec_deg"] = Angle((df[8].astype(float), df[9].astype(float), df[10].astype(float)), unit = u.deg).value
                #df['Mag'] = df[8]
                df['observatory_code'] = observatory_code
                ephemeris_dfs.append(df)

        ephemeris = pd.concat(ephemeris_dfs, ignore_index=True)
        return ephemeris


    def _propagateOrbits(self, orbits, stop, out_dir=None):
        '''
        Propagate orbits to stop time.

        Parameters
        ----------
        orbits : `~orb_it.raiden.Orbits`
            Orbits to propagate.
        stop : `~astropy.time.core.Time`
            Times to which to propagate each orbit.
        out_dir : str, optional
            Save input and output files to this directory. Will create
            a sub directory called propagation inside this directory.

        Returns
        -------
        dfs : `~pandas.DataFrame`
            Orbits propagated to t1
        '''

        with tf.TemporaryDirectory() as tdir:
            prop_dfs=[]
            for i in range(orbits.num_orbits):
                orbit_id_i = orbits.ids[i]
                eq1file(orbits[i], tdir=tdir)
                fopfile(tdir=tdir)
                propHelp(stop[0], tdir=tdir)
        
                home = os.environ["HOME"]
        
                call = [f'{home}/orbfit/src/fitobs/fitobs.x', '<', 'ast.inp']
                ca=' '.join(call)
                with open(os.path.join(tdir,'temp_prop.txt'),'w') as fi:
                    subprocess.call(ca,shell=True,cwd=tdir,stderr=fi,stdout=subprocess.DEVNULL)
                
                # NEW TABLE GEN
                
                ret = open(os.path.join(tdir,'temp_prop.txt')).read().split('\n')
                for j in range(len(ret)):
                    if 'CAR' in ret[j]:
                        p1 = ret[j].split()
                    if 'MJD' in ret[j]:
                        mjd = ret[j].split()
                prop = np.array(p1[1:],dtype=np.float64)
                p_df = pd.DataFrame()
                p_df['orbit_id'] = [orbit_id_i]
                p_df['mjd_tdb'] = Time([mjd[1]],format='mjd',scale='tt').tdb.mjd #returns in tdt
                p_df['x'] = [prop[0]]
                p_df['y'] = [prop[1]]
                p_df['z'] = [prop[2]]
                p_df['vx'] = [prop[3]]
                p_df['vy'] = [prop[4]]
                p_df['vz'] = [prop[5]]
                prop_dfs.append(p_df)
            dfs = pd.concat(prop_dfs, ignore_index=True)   
            return dfs

    def _orbitDetermination(self, observations, out_dir=None):
        '''
        Determine orbits from a set of observations.

        Parameters
        ----------
        observations : `~pandas.DataFrame`
            DataFrame containing the observational data.
        full_output : `bool`, optional
            If True, return the ephemeris and the process return codes and residuals.
        out_dir : `str`, optional
            Path to the output directory.
        
        Returns
        -------
        od_orbits : `~pandas.DataFrame`
            DataFrame containing the orbital data.
        '''
        _observations = observations.copy()
        _observations.sort_values(
                by=["mjd_utc"],
                inplace=True
            )
        id_present = False
        id_col = None
        if "permID" in _observations.columns.values:
            id_present = True
            id_col = "permID"
        if "provID" in _observations.columns.values:
            id_present = True
            id_col = "provID"
        if "trkSub" in _observations.columns.values:
            id_present = True
            id_col = "trkSub"
        if not id_present:
            _observations["trkSub"] = _observations["orbit_id"]
            id_col = "trkSub"
        od_dfs = []
        for i, orbit_id in enumerate(_observations[id_col].unique()):
            orbit_id_short = f"o{i:07d}"
            if "orbit_id" in _observations.columns:
                orbit_id_i = orbit_id
                orbit_id_path = "_".join(orbit_id.split(" "))
                orbit_id_path = "_".join(orbit_id_path.split("/"))
            else:
                orbit_id_i = orbit_id_short
            with tf.TemporaryDirectory() as tdir:
                if 0 not in _observations.columns:
                    print('b')
                    ra = Angle(_observations['RA_deg'],unit='deg').hms
                    dec = Angle(_observations['Dec_deg'],unit='deg').hms
                    day = []
                    mon = []
                    yr = []
                    for i in range(len(_observations)):
                        ti = float(_observations['mjd_utc'][i])
                        mod = ti %1.0
                        sa = str(mod).lstrip('0')
                        date = Time(ti,format="mjd").strftime(f"%d{sa} %b %Y") # NEED TO CHECK THIS FOR INACCURACY
                        day.append(date.split()[0])
                        mon.append(date.split()[1])
                        yr.append(date.split()[2])
                        
                    _observations[1] = mon
                    _observations[0] = day
                    _observations[2] = yr
                    _observations[5] = ra.h
                    _observations[6] = ra.m
                    _observations[7] = ra.s
                    _observations[8] = dec.h
                    _observations[9] = dec.m
                    _observations[10] = dec.s
                os.makedirs(f'{tdir}/obs',exist_ok=True)
                mondic={"Jan":"01",
                "Feb":"02",
                "Mar":"03",
                "Apr":"04",
                "May":"05",
                "Jun":"06",
                "Jul":"07",
                "Aug":"08",
                "Sep":"09",
                "Oct":"10",
                "Nov":"11",
                "Dec":"12"
                }
                # with open(f'{tdir}/obs/temp.obs','w') as f:
                #     for i in range(len(_observations)):
                #         ti = float(_observations['mjd_utc'][i])
                #         mod = ti %1.0
                #         sa = str(mod).lstrip('0')
                #         date = Time(ti,format="mjd").strftime(f"%d{sa}")
                #         day=f"{float(date):.5f}".rjust(8,'0')
                #         if float(_observations[8][i]) > 0:
                #             dc1 = f"+{_observations[8][i]}"
                #         else:
                #             dc1 = _observations[8][i]
                #         st1=f'0temp{"       "}   {_observations[2][i]} {mondic[_observations[1][i]]} {day} {_observations[5][i]} {_observations[6][i]} '+f'{float(_observations[7][i]):.2f}'.rjust(5,'0')+f' {dc1} {_observations[9][i]} '+f'{float(_observations[10][i]):0.2f}'.rjust(5,"0")+f'          {"     "}      {_observations["observatory_code"][i]}\n'
                #         f.write(st1)
                #     f.close()
                band = '    '
                if 11 in _observations.columns:
                    band = 'V   '
                elif 'magV' in _observations.columns:
                    _observations[11] = _observations['magV']
                    band = 'V   '
                else:
                    _observations[11] = "0"
                      
                #hs = '# version=2017\npermID |mode|stn |prog|obsTime                 |ra                  |dec                 |astCat  |mag  |band|photCat |disc|subFmt|precTime|precRA|precDec|notes\n'
                hs = '# version=2017\ntrkSub|obsTime|ra|dec|rmsRA|rmsDec|mag|rmsMag|band|stn|mode|astCat\n'
                with open(f'{tdir}/obs/temp.obs_psv','w') as f:
                    f.write(hs)
                    for i in range(len(_observations)):
                        if float(_observations['Dec_deg'][i]) > 0:
                            dc1 = f"{_observations['Dec_deg'][i]:0.16f} "
                        else:
                            dc1 = f"{_observations['Dec_deg'][i]:0.16f}"
                        if '.' not in  _observations[11][i]:
                            mag = "     "
                        else:
                            mag = f"{_observations[11][i]}".ljust(5,'0')
                        #fs = f"temp   | CCD|{_observations['observatory_code'][i]} |    |{Time(_observations['mjd_utc'][i],format='mjd',precision=4).isot}|{_observations['RA_deg'][i]:0.16f}|{dc1}|        |"+f"{mag}"+f"|{band}|        |    |      |        |      |       |\n"
                        fs = f"temp|{Time(_observations['mjd_utc'][i],format='mjd',precision=4).isot}|{_observations['RA_deg'][i]:0.16f}|{_observations['Dec_deg'][i]:0.16f}|{_observations['RA_sigma_mas'][i]*1000:0.18f}|{_observations['Dec_sigma_mas'][i]*1000:0.18f}|{_observations[11][i]}|0.0|{band}|{_observations['observatory_code'][i]}|CCD|\n"
                        f.write(fs)
                    f.close()
                start = _observations['mjd_utc'][0]
                end = _observations['mjd_utc'][-1:].values[0]
                if end-start > self.arc_limit:
                    os.makedirs(f'{tdir}/obs1',exist_ok=True)
                    with open(f'{tdir}/obs1/temp.obs_psv','w') as f:
                        f.write(hs)
                        i=0
                        while _observations['mjd_utc'][i]-start<self.arc_limit+5 and i<8:
                            if float(_observations['Dec_deg'][i]) > 0:
                                dc1 = f"{_observations['Dec_deg'][i]:0.16f} "
                            else:
                                dc1 = f"{_observations['Dec_deg'][i]:0.16f}"
                            if '.' not in  _observations[11][i]:
                                mag = "     "
                            else:
                                mag = f"{_observations[11][i]}".ljust(5,'0')
                            #fs = f"temp   | CCD|{_observations['observatory_code'][i]} |    |{Time(_observations['mjd_utc'][i],format='mjd',precision=4).isot}|{_observations['RA_deg'][i]:0.16f}|{dc1}|        |"+f"{mag}"+f"|{band}|        |    |      |        |      |       |\n"
                            fs = f"temp|{Time(_observations['mjd_utc'][i],format='mjd',precision=4).isot}|{_observations['RA_deg'][i]:0.16f}|{_observations['Dec_deg'][i]:0.16f}|{_observations['RA_sigma_mas'][i]*1000:0.18f}|{_observations['Dec_sigma_mas'][i]*1000:0.18f}|{_observations[11][i]}|0.0|{band}|{_observations['observatory_code'][i]}|CCD|\n"
                            f.write(fs)
                            i+=1
                        f.close()
                    fopOD2(tdir=tdir)
                else:
                    fopOD(tdir=tdir)
                od1Help(tdir=tdir)
                
                home = os.environ["HOME"]

                
                subprocess.call(f'{home}/orbfit/src/fitobs/fitobs.x < ast.inp',cwd=tdir,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
                
                v1 = open(f'{tdir}/temp.fou').read().split('\n')
                for i in range(len(v1)):
                    if 'KepElem' in v1[i]:
                        lin = v1[i]
                        break
                a1=np.array(lin.split()[1:-2],dtype=np.float64)
                a1=np.append(a1,np.float64(lin.split()[-1].replace(')','')))
                # for i in range(len(v1)):
                #     if 'preliminary orbit elements' in v1[i]:
                #         tim = v1[i].split()[-1]
                #         elem = np.array(v1[i+1].split(),dtype=np.float64)
                #         break
                # a1=np.append(elem,np.float64(tim))
                fopOD(tdir=tdir)
                eqOD(a1,tdir)
                od2Help(tdir)
                with open(f'{tdir}/ls_out.txt','w') as f2:
                    subprocess.call(f'{home}/orbfit/src/fitobs/fitobs.x < ast.inp',cwd=tdir,shell=True,stdout=f2,stderr=subprocess.DEVNULL)
                    f2.close()
                
                v2= open(f'{tdir}/ls_out.txt').read().split('\n')
                for i in range(len(v2)):
                    if 'new elem values' in v2[i]:
                        val=v2[i+1]
                
                vec = np.array(val.split(),dtype=np.float64)
                odd = pd.DataFrame()
                odd['orbit_id'] = [orbit_id_i]
                odd['mjd_utc'] = Time(a1[6],format='mjd',scale='tt').utc.mjd
                odd['mjd_tdb'] = Time(a1[6],format='mjd',scale='tt').tdb.mjd # CHECK TIME SCALE utc or tt
                odd['x'] = vec[0]
                odd['y'] = vec[1]
                odd['z'] = vec[2]
                odd['vx'] = vec[3]
                odd['vy'] = vec[4]
                odd['vz'] = vec[5]
                od_dfs.append(odd)
        od_orbits = pd.concat(od_dfs, ignore_index=True)
        return od_orbits 