from setuptools import setup

setup(
   packages=find_packages(
      where='..',
      include=['orb_it*']),
   use_scm_version={
      "write_to": "orb_it/version.py",
      "write_to_template": "__version__ = '{version}'",
    }
)
