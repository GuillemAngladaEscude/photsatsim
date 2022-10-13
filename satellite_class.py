from orbit_class import Orbit
from attitude_class import Attitude

import astropy.time.core
import numpy as np
from scipy.optimize import root
import scipy.constants as constant
from numpy import pi, ndarray
from astropy.time import Time

from astropy.utils import iers
#iers.conf.iers_auto_url = 'https://datacenter.iers.org/data/9/finals2000A.all'

from astropy import units as u
import astropy.coordinates
from astropy.coordinates import SkyCoord
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body, get_sun, get_body_barycentric, get_body_barycentric_posvel
solar_system_ephemeris.set('jpl') #'de432s'

import aux_functions as auxf

class Satellite:

    def __init__(self, optic, t0, zp: float, za: float, i: float, raan: float, theta_0 = 0, aop = 0):

        self.orbit = Orbit(t0, zp, za, i, raan, theta_0, aop)
        self.att = Attitude(self.orbit)
        self.optic = optic


    def get_stars_in_frame(self,stardata): # aquest mètode va aquí perquè depèn tant d'attitude com d'òptica
        stars_in_frame = []
        for star in stardata:

            # coordenades ICRS de cada estel:
            ra_ICRS = star[2]
            dec_ICRS = star[3]

            # vector que apunta a l'estel en la base Satellite Based Field Of View Fixed:
            r_SBFOVF = self.att.ICRS_2_SBFOVF(dec_ICRS, ra_ICRS)
            r_pqr = self.att.ICRS_2_pqr(dec_ICRS, ra_ICRS)
            u, v = r_pqr[0], r_pqr[1]

            # coordenades de l'estel en píxels dins del frame (new):
            if self.optic.check_if_inframe(r_SBFOVF):
                frame_coord = self.optic.project_SBFOVF_2_sensor(r_SBFOVF, 'linear')
                new_row = [frame_coord[0], frame_coord[1],
                           float(star[5])]  # star[5] dona la magnitud aparent; canviar-ho a star[4] si es vol el flux
                stars_in_frame.append(new_row)

        return stars_in_frame

