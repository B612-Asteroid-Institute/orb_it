import os

__all__ = [
    "Backend"
]

# Will add workers here if multiprocessing is available
# To whoever wants multi-processing support: Look at Joachim Moeyens' THOR backend.py file for inspiration.

class Backend:

    def __init__(self, name="Backend", **kwargs):
        self.__dict__.update(kwargs)
        self.name = name
        self.is_setup = False
        return

    def setup(self):
        return

    def _propagateOrbits(self, orbits, t1):
        """
        Propagate orbits from t0 to t1.

        THIS FUNCTION SHOULD BE DEFINED BY THE USER.

        """
        err = (
            "This backend does not have orbit propagation implemented."
        )
        raise NotImplementedError(err)

    # add mp propagateOrbits()

    def _generateEphemeris(self, orbits, observers):
        """
        Generate ephemerides for the given orbits as observed by
        the observers.

        THIS FUNCTION SHOULD BE DEFINED BY THE USER.

        """
        err = (
            "This backend does not have ephemeris generation implemented."
        )
        raise NotImplementedError(err)

    # add mp generateEphemeris()

    def _orbitDetermination(self):
        err = (
            "This backend does not have orbit determination implemented."
        )
        raise NotImplementedError(err)

    # add mp orbitDetermination()