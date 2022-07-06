import os

__all__ = [
    "Backend"
]

# Will add workers here if multiprocessing is available
# To whoever wants multi-processing support: Look at Joachim Moeyens' THOR backend.py file for inspiration.

class Backend:
    '''
    Basic backend class.

    This was created as a backup and template for new integrators to be made for orb_it or other programs.
    '''
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
        """
        Fit ephemerides observations to determine orbit state.

        THIS FUNCTION SHOULD BE DEFINED BY THE USER.
        
        """
        err = (
            "This backend does not have orbit determination implemented."
        )
        raise NotImplementedError(err)

    # add mp orbitDetermination()