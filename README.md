# Orbital Integrator Tester (Orb_it)

**Orb_it** is a soon-to-be end-to-end testing package that will be able to test orbit prediction accuracy for Open Source Orbit Integrators such as [Find_Orb](https://www.projectpluto.com/find_orb.htm), [OpenOrb](https://github.com/oorb/oorb), and [OrbFit](http://adams.dm.unipi.it/~orbmaint/orbfit/OrbFit/doc/help.html). It currently *only* works for Find_Orb and OpenOrb on Unix based systems. This is an expansion of the earlier testing platform [Validate Find_Orb](https://github.com/B612-Asteroid-Institute/validate_findorb), and shares many features with it.


## Installation and Dependencies

**This project requires** a Linux Based System, Python 3, and working *current* installation of each integrator (preferably installed with conda if possible). It is beneficial to also have a [conda based package manager](https://docs.conda.io/en/latest/) and [Jupyter](https://jupyter.org/) for ease of use.

To install *Orb_it* as a python package, clone this repository, and `cd` into the folder. Then follow the directions below for your use case,

### Installing Dependencies with conda
To install its dependencies in a new conda environment, use this command,
```
conda create -n myenv -c defaults -c conda-forge --file requirements.txt python=3.9
```
---
To install on a preexisting conda environment, first activate your environment, `conda activate myenv`, then,
```
conda install -c defaults -c conda-forge --file requirements.txt
```
---
Once all the dependencies have been installed, run this command to install the development version of this package,
```
python setup.py develop --no-deps
```
### Installing Dependencies with pip
Just type in the command line,
```
python setup.py develop
```
**Note:** It is not recommended to use `python setup.py install` for installation since this project is still in development and may have frequent updates.

---
Then check if its properly installed by typing in the python command line,
```
>>> import orb_it
```

BIG NOTE:
**This project is still in REALLY ACTIVE development** and will have full documentation and a working demo soon!

If you have any more questions email me here: aidan@b612foundation.org
