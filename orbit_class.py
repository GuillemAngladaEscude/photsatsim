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

import aux_functions as auxf

R = 6378.1          # [km] - Earth's equatorial radius (not mean!)
mu = 3.986004418e5  # [km^3/s^2] - Earth's gravitational parameter
J2 = 1.08263e-3     # [] - Earth's second zonal harmonic

class Orbit:
    """
        A class used to define the desired orbital configuration of a satellite at a given instant and propagate it to a
         defined instant delta_t minutes in the future.

        An Orbit can be initialised either with .from_elements() or .from_rv(). The former requires the set of six
         keplerian orbital elements whereas the latter requires the position and velocity vectors at the initial
         instant. Both require the initial time t0 to be introduced as a Time object <class 'astropy.time.core.Time'>

        This class uses the modules numpy and astropy.

        Attributes
        ----------
        t: astropy Time object <class 'astropy.time.core.Time'>
            Current simulation time after the last propagation

        ground_stations: list of str
            Names of the added ground stations

        gs_contact: boolean list
            List of boolean variables containing True if there is visual contact between the satellite and the ground
            station corresponding to the same position in the list ground_stations at the current instant and False
            otherwise

        r_geo_inert: numpy array (float64)
            Position vector of the satellite at the current instant in the non-rotating (inertial) geocentric equatorial
            frame [km, km, km]

        v_geo_inert: numpy array (float64)
            Velocity vector of the satellite at the current instant in the non-rotating (inertial) geocentric equatorial
            frame [km/s, km/s, km/s]

        r_geo_rot: numpy array (float64)
            Position vector of the satellite at the current instant in the rotating (non-inertial) geocentric equatorial
            frame [km, km, km]

        lon: float64
            Longitude of the sub-satellite point at the current instant [deg - range: 0 to 360]

        lat: float64
            Latitude of the sub-satellite point at the current instant [deg - range: 0 to 360]

        altitude: float64
            Altitude of the satellite over the Earth's surface, assuming spherical Earth of R = 6378.1km (equatorial)
            [km]

        Methods
        -------
        from_elements(t0, zp, za, i, raan, theta_0, aop):
            Initialise an Orbit object with th set oif six keplerian elements at the initial instant

        from_rv(t0, r, v):
            Initialise an Orbit object with the inertial geocentric equatorial position and velocity vectors at the
            initial instant.

        add_ground_station(GS_name):
            Adds a Ground Station to the Orbit object with which visual contact will be checked at each propagation

        print_info():
            Prints the orbital information after the last executed propagation.

        propagate(delta_t, method):
            Propagates the position and parameters of the satellite one delta_t into the future with the selected
            method. If no method is provided, "Kepler_J2" is used by default.

        rv_2_elements(r, v) (static method):
            Calculate the six keplerian orbital elements given a position and a velocity vector in the non-rotating
            geocentric equatorial frame.

    """

    all = []

    def __init__(self, t0, zp: float, za: float, i: float, raan: float, theta_0=0, aop=0):
        """
        Initialise an Orbit object
        :param t0: time of the initial conditions of the satellite as a Time object <class 'astropy.time.core.Time'>
        :param zp: perigee altitude above the Earth's surface, considering R = 6378.1km [km]
        :param za: apogee altitude above the Earth's surface, considering R = 6378.1km [km]
        :param i: orbit inclination [deg]
        :param raan: Right Ascension of the Ascending Node [deg]
        :param theta_0: True anomaly [deg]
        :param aop: Argument of perigee [deg]
        """

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
        self.ground_stations = []
        self.gs_contact = []

        # initialise the propagated parameters
        self.propagate(0)
        print("Orbit initialised!")
        self.print_info()

    @classmethod
    def from_elements(cls,t0, zp, za, i, raan, theta_0, aop):
        """
        Initialise an Orbit object with th set of six keplerian elements at the initial instant
        :param t0: time of the initial conditions of the satellite as a Time object <class 'astropy.time.core.Time'>
        :param zp: perigee altitude above the Earth's surface, considering R = 6378.1km [km]
        :param za: apogee altitude above the Earth's surface, considering R = 6378.1km [km]
        :param i: orbit inclination [deg]
        :param raan: Right Ascension of the Ascending Node [deg]
        :param theta_0: True anomaly [deg]
        :param aop: Argument of perigee [deg]
        :return: Initialised Orbit object
        """

        return cls(t0, zp, za, i, raan, theta_0, aop)

    @classmethod
    def from_rv(cls, t0, r, v):
        """
        Initialise an Orbit object with the inertial geocentric equatorial position and velocity vectors at the initial
            instant.
        :param t0: time of the initial conditions of the satellite as a Time object <class 'astropy.time.core.Time'>
        :param r: array with the [x, y, z] components of the satellite's position in the non-rotating Earth-centred,
            equatorial frame [km, km, km]
        :param v: array with the [x, y, z] components of the satellite's velocity in the non-rotating Earth-centred,
            equatorial frame [km/s, km/s, km/s]
        :return: Initialised Orbit object
        """
        [zp, za, i, raan, theta, aop] = cls.rv_2_elements(r, v)
        return cls(t0, zp, za, i, raan, theta, aop)

    def add_ground_station(self,GS_name):
        """
        Adds a Ground Station to the Orbit object with which visual contact will be checked at each propagation
        :param GS_name: String with the name of the desired Ground Station. Current accepted inputs are "Montsec" and
            "Svalbard"
        :return: N/A
        """
        self.ground_stations.append(GS_name)

    def print_info(self):
        """
        Prints the orbital information after the last executed propagation. This information is displayed in 4 blocks:
            - Time: current simulation time, last propagated delta_t, time of last perigee and time since first perigee.
            - Orbital elements: perigee altitude, apogee altitude, semimajor axis, eccentricity, inclination, right
                ascension of the ascending node, argument of perigee, orbital period, true anomaly and specific angular
                momentum.
            - Position: position and velocity vectors in the geocentric equatorial non-rotating frame, position vector
                in the geocentric equatorial rotating frame, longitude, latitude and altitude.
            - Contact with Ground stations: list of the added ground stations, together with a 'Yes' or 'No' depending
                on whether or not the satellite has visual contact with them (assuming 0 degrees over the horizon enough
                to make visual contact)
        :return: N/A
        """

        if self.__e == 0:
            aop_print = 'N/A (circular orbit)'
            tsp_print = 'N/A (circular orbit)'
        else:
            aop_print = f"{self.__w * 180 / pi} degrees"
            tsp_print = f"{(self.__t_0 + self.__propagated_deltat) / 60} minutes"

        if len(self.ground_stations) == 0:
            gs_print = "N/A. Add ground stations using .add_ground_station()"
        else:
            gs_print = ""
            i = 0
            for gs_name in self.ground_stations:
                gs_print += gs_name
                gs_print += ": "
                if self.gs_contact[i]==True:
                    gs_print += "Yes"
                else:
                    gs_print += "No"
                gs_print += "\n    "
                i = i+1

        print(f"""
----------------------------------------------------------------------------------------------------------------------
Current time: {self.t}
Propagated delta_t = {self.__propagated_deltat / 60} minutes:
Time since perigee of initial orbit: {tsp_print};
Time of last perigee: {self.t_last_perigee}
        
Orbital elements at current time:
    Perigee altitude: {self.__a * (1 - self.__e) - R} km
    Apogee altitude: {self.__a * (1 + self.__e) - R} km
    Semimajor axis: {self.__a} km;
    Eccentricity: {self.__e};
    Inclination: {self.__i * 180 / pi} degrees;
    Right Ascension of the Ascending Node: {self.__RAAN * 180 / pi} degrees;
    Argument of perigee: {aop_print};
    Orbital period: {self.__T / 60} minutes;
    True anomaly: {self.__theta * 180 / pi} degrees;
    Specific angular momentum: {self.__h0} km^2/s;
          
Geocentric equatorial inertial frame position: {self.r_geo_inert} km
Geocentric equatorial inertial frame velocity: {self.v_geo_inert} km/s
Geocentric equatorial Earth-fixed frame position: {self.r_geo_rot} km
Longitude: {self.lon} degrees
Latitude: {self.lat} degrees
Altitude: {self.altitude} km

Visual contact with set ground stations:
    {gs_print}
        """)

    @property
    def t_last_perigee(self):
        """
        Returns the day and time at which the satellite passed through the perigee of its orbit for the last time
        :return: time of last perigee as a Time object <class 'astropy.time.core.Time'>
        """
        return self.t - astropy.time.core.TimeDelta(self.__t_0 * u.s) - astropy.time.core.TimeDelta((self.__propagated_deltat - floor(self.__propagated_deltat/self.__T)*self.__T) * u.s)

    @property
    def v_geo_gcrs(self):
        """
        Satellite velocity in the Geocentric Celestial Reference Frame (GCRS)
        :return: Satellite velocity in the Geocentric Celestial Reference Frame (GCRS) [km/s]. At the current iteration,
        assumed to be equal to the velocity in the non-rotating geocentric equatorial frame calculated by
        __get_r_geo_inertial.
        """
        return self.v_geo_inert

    def propagate(self, delta_t, method="Kepler_J2"):
        """
        Propagates the position and parameters of the satellite one delta_t into the future with the selected method. If
            no method is provided, "Kepler_J2" is used by default.
        :param delta_t: desired propagation time interval [min]
        :param method: desired propagation method. Current accepted inputs are:
            - "Kepler_J2": applies a keplerian two-body propagation adding the average secular effect of the J2 zonal
                harmonic on the right ascension of the ascending node and the argument of perigee.
        :return: N/A
        """
        self.t = self.__t0 + astropy.time.core.TimeDelta(delta_t * u.min)
        delta_t = delta_t * 60
        self.__propagated_deltat = delta_t
        if method == "Kepler_J2":
            self.__propagate_KeplerJ2(delta_t)
        elif method == "interpolate_from_file":
            pass  # quan s'hagi propagat una òrbita amb precisió i se n'hagin guardat els valors en un arxiu separat
                  # s'utilitzarà aquest mètode
        else:
            raise Exception(f"Method {method} is not recognized as a valid method!")

        self.__get_lat_lon_h()

    def __propagate_KeplerJ2(self, delta_t):
        """
        Update the propagated orbital elements after delta_t seconds with the Kepler+J2 method (keplerian two-body
            motion adding the average secular effect of the J2 zonal harmonic on the right ascension of the ascending
            node and the argument of perigee).
        :param delta_t: desired propagation time interval [s]
        :return: N/A
        """

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

    def __solveKepler(self, Me):
        """
        Solves Kepler's equation (M_e = E - e*sin(E)) for the eccentric anomaly E given a mean anomaly Me, by means of
            a Newton-Raphson iteration
        :param Me: mean anomaly [rad]
        :return: eccentric anomaly E [rad]
        """

        # initial guess:
        if Me < pi:
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
        """
        Calculates and stores the non-rotating geocentric equatorial position and velocity vectors and their modules
            [km], [km/s]
        :return: N/A (values stored in the class attributes)
        """

        self.__mod_r = self.__h ** 2 / mu * 1 / (1 + self.__e * np.cos(self.__theta))
        r_perif = np.multiply(self.__mod_r, [np.cos(self.__theta), np.sin(self.__theta), 0])
        self.r_geo_inert = np.matmul(self.__Q(), r_perif)

        self.__mod_v = mu / self.__h
        v_perif = np.multiply(self.__mod_v, [-np.sin(self.__theta), self.__e + np.cos(self.__theta), 0])
        self.v_geo_inert = np.matmul(self.__Q(), v_perif)

    def __get_r_geo_rot(self):
        """
        Calculates and stores the Earth-fixed geocentric equatorial position and velocity vectors and their modules
            [km], [km/s]. The angle between the Greenwich Meridian and the Vernal Equinox is obtained by astropy as
            the Greenwich Mean Sidereal Time (GMST) at the current instant.
        :return: N/A (values stored in the class attributes)
        """

        self.__get_r_geo_inertial()
        t0_GMST = self.t.sidereal_time('mean', 'greenwich') #'apparent' instead of 'mean' also accounts for Earth's
                                                            # nutation
        self.r_geo_rot = np.matmul(np.transpose(auxf.R3(t0_GMST.rad)), self.r_geo_inert)
        self.__check_GS_contact()

    def __get_lat_lon_h(self):
        """
        Calculates the (lon,lat) [deg, deg] coordinates of the Sub Satellite Point and its altitude over the surface
            [km], assuming a perfectly spherical Earth of radius R = 6378.1 km (equatorial)
        :return: N/A (values stored in the class attributes)
        """

        self.__get_r_geo_rot()
        l = self.r_geo_rot[0] / self.__mod_r
        m = self.r_geo_rot[1] / self.__mod_r
        n = self.r_geo_rot[2] / self.__mod_r

        lat = np.arcsin(n)
        lon = np.arccos(l / np.cos(lat))
        if m <= 0:
            lon = 2 * pi - lon

        self.lat = lat * 180 / pi
        self.lon = lon * 180 / pi
        self.altitude = self.__mod_r - R

    def __check_GS_contact(self):
        """
        Check if there is visual contact between the satellite at its current position and the defined Ground Stations,
            assuming spherical Earth of R = 6378.1 km and the visual contact limit to be at 0 deg over the horizon.
        :return: N/A (values stored in the attribute gs_contact)
        """

        self.gs_contact = []
        for gs_name in self.ground_stations:
            if gs_name =="Montsec":
                x_GS = 4.729730e+03
                y_GS = 6.023677e+01
                z_GS = 4.266734e+03
                r_GS = np.array([x_GS, y_GS, z_GS])

            elif gs_name == "Svalbard":
                x_GS = 1.252730e+03
                y_GS = 3.452427e+02
                z_GS = 6.236222e+03
                r_GS = np.array([x_GS, y_GS, z_GS])

            else:
                raise Exception(f"Ground station name {gs_name} not acknowledged!")
                # print(f"Ground station name {gs_name} not acknowledged!")
                # exit()

            r_sat = self.r_geo_rot
            d = np.linalg.norm(r_sat - r_GS)
            h = self.altitude
            d_max = np.sqrt(h**2 + 2*h*R)

            if d<d_max:
                self.gs_contact.append(True)
            else:
                self.gs_contact.append(False)


    @staticmethod
    def rv_2_elements(r, v):
        """
        Calculate the six keplerian orbital elements given a position and a velocity vector in the
            non-rotating geocentric equatorial frame (static method).
        :param r: array with the [x, y, z] components of the satellite's position in the non-rotating Earth-centred,
            equatorial frame [km, km, km]
        :param v: array with the [x, y, z] components of the satellite's velocity in the non-rotating Earth-centred,
            equatorial frame [km/s, km/s, km/s]
        :return: array with six orbital elements in the following order: perigee altitude, apogee altitude, inclination,
            right ascension of the ascending node, true anomaly, argument of perigee [km, km, deg, deg, deg, deg]
        """

        r_mod = np.linalg.norm(r)
        v_mod = np.linalg.norm(v)
        a = 1/(2/r_mod - v_mod**2/mu)
        vr = np.dot(r,v)/r_mod
        e_vec = (1/mu)*((v_mod**2 - mu/r_mod)*r - r_mod*vr*v)
        e = np.linalg.norm(e_vec)
        rp = a*(1-e)
        ra = a*(1+e)
        zp = rp-R
        za = ra-R
        h = np.cross(r,v)
        i = np.arccos(h[2]/np.linalg.norm(h))
        N = np.cross([0,0,1],h)
        N_mod = np.linalg.norm(N)
        RAAN = np.arccos(N[0]/N_mod)
        if N[1]<0:
            RAAN = 2*pi-RAAN
        aop = np.arccos((np.dot(N,e_vec))/(N_mod*e))
        if e_vec[2]<0:
            aop = 2*pi-aop
        theta = np.arccos((np.dot(e_vec,r))/(e*r_mod))
        if vr<0:
            theta = -theta
        return np.array([zp, za, i*180/pi, RAAN*180/pi, theta*180/pi, aop*180/pi])

    def __Q(self):
        """
        Calculates the change of basis matrix Q_xX, perifocal to inertial geocentric equatorial, based on the current
            values of inclination i, right ascension of the ascending node (RAAN) and argument of perigee (w)
        :return: change of basis matrix Q_xX, perifocal to inertial geocentric equatorial
        """
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
