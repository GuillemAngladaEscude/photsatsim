import astropy.time.core
import numpy as np
from scipy.optimize import root
import scipy.constants as constant
from numpy import pi, ndarray
from astropy.time import Time

from math import floor

from astropy.utils import iers
# iers.conf.iers_auto_url = 'https://datacenter.iers.org/data/9/finals2000A.all'

from astropy import units as u
import astropy.coordinates
from astropy.coordinates import SkyCoord
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body, get_sun, get_body_barycentric, get_body_barycentric_posvel
solar_system_ephemeris.set('jpl') #'de432s'


class Orbit:
    R = 6378.1
    mu = 3.986004418e5  # [km^3/s^2]
    J2 = 1.08263e-3  # []
    all = []

    def __init__(self, t0, zp: float, za: float, i: float, raan: float, theta_0=0, aop=0):

        R = 6378.0
        mu = 3.986004418e5  # [km^3/s^2]
        J2 = 1.08263e-3  # []
        a = R + (zp + za) / 2
        e = (za - zp) / (2 * R + zp + za)

        # run validations to arguments
        assert zp >= 100, f"Perigee altitude cannot be smaller than 100km!"
        assert za >= zp, f"Apogee altitude cannot be smaller than perigee altitude!"
        assert i <= 180, f"Inclination cannot be greater than 180 degrees!"
        assert i >= 0, f"Inclination cannot be negative!"

        # assign
        self.__t0 = t0
        self.__a0 = a
        self.__i0 = i * pi / 180
        self.__e0 = e
        self.__RAAN0 = raan * pi / 180
        self.__w0 = aop * pi / 180
        self.__theta0 = theta_0 * pi / 180

        self.__RAAN_dot = -(3 / 2 * (mu ** 0.5) * J2 * (R ** 2) / (((1 - e ** 2) ** 2) * (a ** (7 / 2)))) * np.cos(
            self.__i0)
        self.__w_dot = -(3 / 2 * (mu ** 0.5) * J2 * (R ** 2) / (((1 - e ** 2) ** 2) * (a ** (7 / 2)))) * (
                    2.5 * (np.sin(self.__i0)) ** 2 - 2)
        self.__T0 = 2 * pi * a ** 1.5 / np.sqrt(mu)
        self.__h0 = np.sqrt((R + zp) * mu * (1 + self.__e0))

        self.__propagated_deltat = 0.0

        # Eccentric anomaly
        E = 2*np.arctan(np.sqrt((1 - self.__e0) / (1 + self.__e0)) * np.tan(self.__theta0 / 2))
        # Mean anomaly:
        Me = E - self.__e0 * np.sin(E)
        # Time since perigee:
        self.__t_0 = Me * self.__T0 / (2 * pi)

        # actions to execute
        Orbit.all.append(self)  # add to list

        # initialize the propagated parameters
        self.propagate(0)
        self.get_lat_lon_h()
        print("Orbit initialised!")
        self.print_elements()

    def print_elements(self):
        R = 6378.0

        if self.__e == 0:
            aop_print = 'N/A'
            tsp_print = 'N/A'
        else:
            aop_print = f"{self.__w * 180 / pi} degrees"
            tsp_print = f"{(self.__t_0 + self.__propagated_deltat) / 60} minutes"

        print(f"""
        Orbital elements after propagating t = {self.__propagated_deltat / 60} minutes:
          Perigee altitude: {self.__a * (1 - self.__e) - R} km
          Apogee altitude: {self.__a * (1 + self.__e) - R} km
          Semimajor axis: {self.__a} km;
          Eccentricity: {self.__e};
          Inclination: {self.__i * 180 / pi} degrees;
          Right Ascension of the Ascending Node: {self.__RAAN * 180 / pi} degrees;
          Argument of perigee: {aop_print};
          Orbital period: {self.__T / 60} minutes;
          True anomaly: {self.__theta * 180 / pi} degrees;
          Time since perigee: {tsp_print};
          Specific angular momentum: {self.__h0} km^2/s.
          Geocentric equatorial inertial frame position: {self.__r_geo_inert} km
          Geocentric equatorial Earth-fixed frame position: {self.__r_geo_rot} km
          Longitude: {self.__lon} degrees
          Latitude: {self.__lat} degrees
          Altitude: {self.__altitude} km
        """)

    @property
    def t(self):
        return self.__t

    @property
    def t_last_perigee(self):

        return self.__t - astropy.time.core.TimeDelta(self.__t_0 * u.s) #- astropy.time.core.TimeDelta(floor(self.__propagated_deltat/self.__T0)*self.__T0 * u.s)

    def propagate(self, delta_t, method="Kepler_J2"):
        self.__t = self.__t0 + astropy.time.core.TimeDelta(delta_t * u.min)
        #self.__delta_t = delta_t * 60
        self.__propagated_deltat = delta_t * 60
        if method == "Kepler_J2":  # Aquest mètode aplica una simple propagació kepleriana (2-body) afegint-hi l'efecte secular de J2 sobre RAAN i AOP.
            self.__propagate_KeplerJ2(delta_t)
        else:
            print(f"Method {method} is not recognized as a valid method!")

        self.__get_lat_lon_h()

    def __propagate_KeplerJ2(self, delta_t):
        mu = 3.986004418e5

        t = self.__t_0 + delta_t

        # New mean anomaly:
        Me = 2 * pi * t / self.__T0
        # New eccentric anomaly: Solve Kepler equation
        E = self.__solveKepler(Me)

        self.__theta = 2 * np.arctan(np.sqrt((1 + self.__e0) / (1 - self.__e0)) * np.tan(E / 2))
        self.__RAAN = self.__RAAN0 + delta_t * self.__RAAN_dot
        self.__w = self.__w0 + delta_t * self.__w_dot
        self.__a = self.__a0
        self.__e = self.__e0
        self.__i = self.__i0
        self.__h = np.sqrt((self.__a * (1 - self.__e)) * mu * (1 + self.__e))
        self.__T = 2 * pi * self.__a ** 1.5 / np.sqrt(mu)
        #self.__t_0 = delta_t - np.floor(delta_t / self.__T0) * self.__T0  # re-definin temps inicial al temps transcorregut des de l'útlim perigeu

    def __solveKepler(self, Me):

        # initial guess:
        if (Me < pi):
            E0 = Me + self.__e0 / 2
        elif Me > pi:
            E0 = Me - self.__e0 / 2
        elif Me == pi:
            E0 = Me

        delta = 1e-8
        ratio = 1
        E = E0
        while ratio > delta:
            f = E - self.__e0 * np.sin(E) - Me
            fp = 1 - self.__e0 * np.cos(E)
            ratio = f / fp
            E = E - ratio
        return E

    def __get_r_geo_inertial(self):  # geocentric equatorial inertial frame
        mu = 3.986004418e5
        self.__mod_r = self.__h ** 2 / mu * 1 / (1 + self.__e * np.cos(self.__theta))
        r_perif = np.multiply(self.__mod_r, [np.cos(self.__theta), np.sin(self.__theta), 0])
        self.__r_geo_inert = np.matmul(self.__Q(), r_perif)

        self.__mod_v = mu / self.__h
        v_perif = np.multiply(self.__mod_v, [-np.sin(self.__theta), self.__e + np.cos(self.__theta), 0])
        self.__v_geo_inert = np.matmul(self.__Q(), v_perif)

    def __get_r_geo_rot(self):
        wE = 2 * pi / (23.934469 * 3600)
        self.__get_r_geo_inertial()
        theta = wE * self.__propagated_deltat # todo: que depengui de l'hora inicial
        r_geo_rot = np.matmul(np.transpose(self.__R3(theta)), self.__r_geo_inert)
        self.__r_geo_rot = r_geo_rot
        # return r_geo_rot

    def __get_lat_lon_h(self):
        R = 6378.0  # Radi equatorial, no mitjà!
        self.__get_r_geo_rot()
        l = self.__r_geo_rot[0] / self.__mod_r
        m = self.__r_geo_rot[1] / self.__mod_r
        n = self.__r_geo_rot[2] / self.__mod_r

        lat = np.arcsin(n)
        lon = np.arccos(l / np.cos(lat))
        if m <= 0:
            lon = 2 * pi - lon

        self.__lat = lat * 180 / pi
        self.__lon = lon * 180 / pi
        self.__altitude = self.__mod_r - R
        # return [self.__lat, self.__lon, self.__altitude]

    def get_r_geo_inertial(self):
        return self.__r_geo_inert

    def get_r_geo_rot(self):
        return self.__r_geo_rot

    def get_lat_lon_h(self):
        return [self.__lat, self.__lon, self.__altitude]

    def __R3(self, theta):  # todo: revisar les tres matrius. Al curtis surten les transposades com si no ho fossin
        R = np.zeros((3, 3))
        R[0, 0] = np.cos(theta)
        R[0, 1] = -np.sin(theta)
        R[1, 0] = np.sin(theta)
        R[1, 1] = np.cos(theta)
        R[2, 2] = 1
        return R

    def __R2(self, theta):
        R = np.zeros((3, 3))
        R[0, 0] = np.cos(theta)
        R[0, 2] = -np.sin(theta)
        R[1, 1] = 1
        R[2, 0] = np.sin(theta)
        R[2, 2] = np.cos(theta)
        return R

    def __R1(self, theta):
        R = np.zeros((3, 3))
        R[0, 0] = 1
        R[1, 1] = np.cos(theta)
        R[1, 2] = np.sin(theta)
        R[2, 1] = -np.sin(theta)
        R[2, 2] = np.cos(theta)
        return R

    def Rs(self, s, alpha):
        s = s / np.linalg.norm(s)
        sx = s[0]
        sy = s[1]
        sz = s[2]
        Q = np.zeros((3, 3))
        Q[0, 0] = np.cos(alpha) + (1 - np.cos(alpha)) * sx ** 2
        Q[0, 1] = (1 - np.cos(alpha)) * sx * sy - sz * np.sin(alpha)
        Q[0, 2] = (1 - np.cos(alpha)) * sx * sz + sy * np.sin(alpha)
        Q[1, 0] = (1 - np.cos(alpha)) * sx * sy + sz * np.sin(alpha)
        Q[1, 1] = np.cos(alpha) + (1 - np.cos(alpha)) * sy ** 2
        Q[1, 2] = (1 - np.cos(alpha)) * sy * sz - sx * np.sin(alpha)
        Q[2, 0] = (1 - np.cos(alpha)) * sx * sz - sy * np.sin(alpha)
        Q[2, 1] = (1 - np.cos(alpha)) * sy * sz + sx * np.sin(alpha)
        Q[2, 2] = np.cos(alpha) + (1 - np.cos(alpha)) * sz ** 2

    def __Q(self):
        i = self.__i
        RAAN = self.__RAAN
        w = self.__w
        Q = np.zeros((3, 3))
        Q[0, 0] = -np.sin(RAAN) * np.cos(i) * np.sin(w) + np.cos(RAAN) * np.cos(w)
        Q[0, 1] = -np.sin(RAAN) * np.cos(i) * np.cos(w) - np.cos(RAAN) * np.sin(w)
        Q[0, 2] = np.sin(RAAN) * np.sin(i)
        Q[1, 0] = np.cos(RAAN) * np.cos(i) * np.sin(w) + np.sin(RAAN) * np.cos(w)
        Q[1, 1] = np.cos(RAAN) * np.cos(i) * np.cos(w) - np.sin(RAAN) * np.sin(w)
        Q[1, 2] = -np.cos(RAAN) * np.sin(i)
        Q[2, 0] = np.sin(i) * np.sin(w)
        Q[2, 1] = np.sin(i) * np.cos(w)
        Q[2, 2] = np.cos(i)
        return Q
