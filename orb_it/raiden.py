import numpy as np
import pandas as pd
from astropy.time import Time
from astropy import units as u
from astroquery.jplhorizons import Horizons

CARTESIAN_COLS = ["x", "y", "z", "vx", "vy", "vz"]
CARTESIAN_UNITS = [u.au, u.au, u.au, u.au / u.d,  u.au / u.d,  u.au / u.d]
REQ_COLS = CARTESIAN_COLS

class Orbits():
    def __init__(self, data,
                 epochs=None,
                 ids=None,
                 H=None,
                 G=None):
        if type(data) is str:
            self.data = pd.read_csv(data,float_precision="round_trip")
        elif type(data) is not pd.DataFrame:
                raise TypeError("Input must be a file path string to a csv or pandas DataFrame")
        else:
            self.data = data
      
        if np.all(pd.Series(REQ_COLS).isin(self.data.columns)):
            pass
        else:
            raise ValueError("Input data must contain the following columns: {}".format(REQ_COLS))
        if 'epoch' in self.data.columns:
            self.epochs = Time(self.data['epoch'], scale='tdb', format='mjd')
        elif 'mjd_tdb' in self.data.columns:
            self.epochs = Time(self.data['mjd_tdb'], scale='tdb', format='mjd')
        elif epochs is not None:
            self.epochs = epochs
        else:
            raise ValueError("Input data must contain either 'epoch' or 'mjd_tdb' time column")
        
        self.num_orbits = self.data.shape[0]
        
        if 'orbit_id' in self.data.columns:
            self.ids = self.data['orbit_id'].to_numpy(dtype='<U60')
        elif ids is not None:
            self.ids = ids
        else:
            self.ids = f'{np.arange(self.num_orbits)}'
        if 'H' in self.data.columns:
            self.H = self.data['H'].values
        elif H is not None:
            self.H = H
        else:
            self.H = None
        if 'G' in self.data.columns:
            self.G = self.data['G'].values
        elif G is not None:
            self.G = G
        else:
            self.G = None
        self.cartesian = self.data[CARTESIAN_COLS].values
        self.cartesian_units = [u.au, u.au, u.au, u.au / u.d, u.au / u.d, u.au / u.d]
        df = pd.DataFrame({'orbit_id': self.ids,'mjd_tdb': self.epochs.value})
        df[CARTESIAN_COLS] = self.cartesian
        df["H"] = self.H
        df["G"] = self.G
        self.df = df
        return

    def __len__(self):
        return self.num_orbits

    def __getitem__(self, i):
        if isinstance(i, slice):
            return Orbits(self.df[i])
            #if abs(i.start) >= self.num_orbits:
            #    raise IndexError
            #if abs(i.stop) > self.num_orbits:
            #    raise IndexError
        
        elif isinstance(i, int):
            if i < 0:
                i += self.num_orbits
            try:
                self.df.loc[i]
            except:
                raise IndexError

            i = slice(i, i+1)
            if abs(i.start) >= self.num_orbits and abs(i.stop) > self.num_orbits:
                raise IndexError

        return Orbits(self.df[i])

def loadOrb(data):
    '''
    Load orbit data from a .csv file or pandas DataFrame.

    Parameters
    ----------
    data: str or pandas.DataFrame
        Path to the .csv file or pandas DataFrame containing the orbit data.
        Must have the following columns:
        x: float
            x element of the state vector in Astronomical Units.
        y: float
            y element of the state vector in Astronomical Units.
        z: float
            z element of the state vector in Astronomical Units.
        vx: float
            x velocity element of the state vector in Astronomical Units per day.
        vy: float
            y velocity element of the state vector in Astronomical Units per day.
        vz: float
            z velocity element of the state vector in Astronomical Units per day.
        epoch or mjd_tdb: float, also can be a astropy.time.core.Time object converted to float
            Time of the state vector in mjd with a tdb scale.
    
    Returns
    -------
    orbits : orbit object `~validate_findorb.raiden.Orbits`
    '''
    return Orbits(data)

def getOrbHorizons(target, t0):
    '''
    Gets the orbital state vector from JPL Horizons for a given target and time.

    Parameters
    ----------
    target : str or list of str
        Name(s) of the target to get the orbital state vector for.
    t0 : astropy.time.core.Time
        Time object with scale='tdb' format='mjd' for the time of the state vector.

    Returns
    -------
    orbits: orbit object `~validate_findorb.raiden.Orbits`
    '''
    
    if type(target) == str:
        target = [target]
    targets_i = []
    for i in target:
        hobj = Horizons(id=i,epochs=t0.tdb.mjd,location='@sun').vectors(refplane="ecliptic",aberrations="geometric",)
        targets_i.append(hobj.to_pandas())
    targets = pd.concat(targets_i, ignore_index=True)
    return Orbits(targets,ids=targets['targetname'].values,epochs=t0+np.zeros(len(target)))