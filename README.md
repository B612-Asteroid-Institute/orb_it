# Comparing Open Source Minor-planet Orbits (COSMO)

**Cosmo** is a soon-to-be end-to-end testing package that will be able to compare and test orbit prediction accuracy for Open Source Orbit Integrators such as Find_Orb. It currently *only* works for Find_Orb on Unix based systems.


## Installation and Dependencies

### To Install
Clone this repository,

```
git clone https://github.com/berres2002/cosmo.git
```

[THOR]('https://github.com/moeyensj/thor') is currently used for some of its utilities, so that needs to also be downloaded:
```git clone git@github.com:moeyensj/thor.git``` 

Or, 
```git clone https://github.com/moeyensj/thor.git``` 

Then install THOR via:
```
conda activate adam_validation_py38
cd {DIRECTORY WHERE THOR DOWNLOADED}
conda install -c defaults -c conda-forge -c moeyensj --file requirements.txt
python setup.py develop --no-deps
```
