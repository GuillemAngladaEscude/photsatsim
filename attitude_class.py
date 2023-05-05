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
from astropy.coordinates import solar_system_ephemeris, EarthLocation, CartesianRepresentation
from astropy.coordinates import get_body, get_sun, get_body_barycentric, get_body_barycentric_posvel
solar_system_ephemeris.set('jpl') #'de432s'

from orbit_class import Orbit
from optic_class import Optic

import aux_functions as auxf

from warnings import warn

c = 299792.458  # [km/s]
R = 6378.1  # [km]

class Attitude:

    """
    A class used to set, calculate and store the attitude configuration of the satellite and siderostat pointing angle
    at the current epoch. The desired configuration can be set in three ways: 1) specifying a target location in the
    ICRS frame, which calculates the necessary angles to put it in the center of the frame with the current satellite
    orbital state; 2) specifying the desired attitude and siderostat angles, or 3) specifying the desired pointing
    vector in the Satellite-Based, Rotating, Ecliptic-Fixed (SBREF) frame. Methods 2 and 3 compute the resulting target
    in the ICRS frame, and all three methods check for possible Earth, Sun or solar panel accidental screenings.

    This class uses the modules numpy, astropy and warnings.

    Attributes (read-only, accessible as properties)
    ----------
    r_SS_SBREF: double array
        Unit vector indicating the direction of the central pixel in the field of view, in the Satellite-
        Based, Rotating, Ecliptic-Fixed frame (SBREF) [ , , ]

    angles: double array
        Array with the current three attitude angles, yaw, pitch and roll, in this order [deg, deg, deg]

    yaw: double
        Current Yaw angle [deg]

    pitch: double
        Current Pitch angle [deg]

    roll: double
        Current Roll angle [deg]

    r_SS_GCRS: double array
        Unit vector in the Geocentric Celestial Reference Frame (GCRS) frame pointing to the sky coordinates of the
        central pixel of the sensor [ , , ]

    dec_ICRS: double
        ICRS declination of the sky point seen in the central pixel of the frame with the current attitude
        and orbital configuration [deg]

    ra_ICRS: ICRS right ascension of the sky point seen in the central pixel of the frame with the current attitude
        and orbital configuration [deg]

    Methods
    -------
    set_r_SS_SBREF(vec)
            Sets the attitude (yaw, pitch, roll) and siderostat pointing (alpha) parameters to make the sky coordinates
            expressed by a unit vector in the Satellite-Based, Rotating, Ecliptic-Fixed (SBREF) frame lay in the center
            of the frame, and warns of Earth, Sun or solar panels screening if vision blocking is found.

    set_angles(yaw, pitch, roll, alpha)
            Sets the attitude parametres of the satellite to the given values of yaw, pitch and roll angles, and the
            siderostat position to the given alpha angle. Calculates the sky coordinates to which the siderostat points
            under these conditions (accessible through dec_ICRS() and ra_ICRS() properties) and  warns of Earth, Sun or
            solar panels screening if vision blocking is found.

    set_pointing(ra_ICRS, dec_ICRS)
            Sets the attitude (yaw, pitch, roll) and siderostat pointing (alpha) parameters to make the sky coordinates
            expressed by a right ascension and a declination in the International Celestial Reference System (ICRS) lay
            in the center of the frame, and warns of Earth, Sun or solar panels screening if vision blocking is found.

    ICRS_2_SBFOVF(ra_ICRS, dec_ICRS)
            Converts the sky coordinates given by a right ascension and a declination in the International Celestial
            Reference System (ICRS) to a unit vector in the Satellite-Based, Field-Of-View-Fixed (SBFOVF) frame, based
            on the satellite's current orbital and attitude configuration

    ICRS_2_pqr(ra_ICRS, dec_ICRS)
            Converts the sky coordinates given by a right ascension and a declination in the International Celestial
            Reference System (ICRS) to a unit vector in the pqr frame, based on the satellite's current orbital and
            attitude configuration

    SBFOVF_2_pqr(r_SBFOVF)
            Converts the sky coordinates given by a unit vector in the Satellite-Based, Field-of-View-Fixed (SBFOVF) frame
            to a unit vector in the pqr frame. Since both reference frames are attached to the optical field of view,
            this transformation is independent of the orbital and attitude parameters.

    GCRS_2_SBFOVF(r_GCRS)
            Converts the sky coordinates given by a cartesian unit vector in the Geocentric Celestial Reference System
            (GCRS) to its coordinates in the Satellite-Based, Field-Of-View-Fixed (SBFOVF) frame

    get_sector()
            Gets the folder name with the Gaia star data of the sky sector where the siderostat is currently pointing. At
            the current version it always outputs 'test_sector_2'.

    """

    def __init__(self, orbit: Orbit, optic: Optic):

        """
        Initialise the class Attitude. The attitude is by default set to the parking configuration, i.e., yaw=0,
            pitch=0, roll=0, alpha = 180 deg.
        :param Orbit orbit: Orbit object with the desired orbital parameters for the satellite
        :param Optic optic: Optic object with the desired optical parameters for the instrument on board
        """

        self.__solar_panel_screening_angle = 90*pi/180  # hard-coded, TBD

        self.__optic = optic
        self.__orbit = orbit
        self.__t = orbit.t

        #parking values
        yaw_0 = 0
        pitch_0 = 0
        roll_0 = 0
        alpha_0 = 180
        self.set_angles(yaw_0, pitch_0, roll_0, alpha_0)


        self.__angles_2_r_SS_SBREF()
        self.__r_SS_SBREF_2_GCRS()
        self.__r_SS_GCRS_2_angles_ICRS()

    def __get_M_SBREF2GCRS(self):
        """
        Calculates the change-of-basis matrix between the Satellite-Based, Rotating, Ecliptic-Fixed (SBREF) frame and
            the Geocentric Equatorial Reference System (GCRS), based on the orbital parameters at the current epoch.
        :return: Matrix M to change SBREF to GCRS (i.e. the three SBREF column vectors in the GCRS frame), based on the
            orbital parameters at the current epoch. [[ , , , ],[ , , , ],[ , , , ]]
        """

        t = self.__orbit.t_last_perigee
        M = np.zeros((3, 3))

        EM_barycenter_GCRS = get_body('earth-moon-barycenter', t)  # GCRS: Geocentric Celestial Reference Frame
        EM_barycenter_GCRS.representation_type = 'cartesian'
        r = [EM_barycenter_GCRS.x.value, EM_barycenter_GCRS.y.value, EM_barycenter_GCRS.z.value]
        loc = EarthLocation.from_geocentric(r[0], r[1], r[2], u.km)

        k_ecl_GCRS = self.__get_k_ecl_GCRS(t, loc)

        M[0,2] = k_ecl_GCRS[0]
        M[1,2] = k_ecl_GCRS[1]
        M[2,2] = k_ecl_GCRS[2]

        r_sun_unit = self.__get_r_sun_unit_GCRS(t, loc)
        M[0,1] = r_sun_unit[0]
        M[1,1] = r_sun_unit[1]
        M[2,1] = r_sun_unit[2]

        x = np.cross(r_sun_unit, k_ecl_GCRS)
        M[0,0] = x[0]
        M[1,0] = x[1]
        M[2,0] = x[2]

        #print(f"SBREF base determinant: {np.linalg.det(M)}")

        return M

    def __get_k_ecl_GCRS(self, t, loc):
        """
        Outputs the cartesian coordinates of a unit vector pointing to the ecliptic pole in the GCRS frame given an
            instant t and a location loc # todo: confirmar si independent de t i loc i if so treure-ho
        :param t: current epoch as a Time object <class 'astropy.time.core.Time'>
        :param loc: desired origin location as an Earth Location object <class 'astropy.coordinates.earth.EarthLocation'>
        :return: Unit vector pointing to the ecliptic pole in the GCRS frame [ , , , ]
        """

        a = EarthLocation.get_gcrs_posvel(loc, t)
        r_EMBC_GCRS = a[0]  # EMBC: Earth-Moon Barycenter
        v_EMBC_GCRS = a[1]

        BME_frame = astropy.coordinates.builtin_frames.BarycentricMeanEcliptic()

        a = SkyCoord(lon=0, lat=90, obstime=t, frame=BME_frame, unit='deg')
        b = a.transform_to('gcrs')

        c = astropy.coordinates.GCRS(ra=b.ra, dec=b.dec, obstime=t, obsgeoloc=r_EMBC_GCRS, obsgeovel=v_EMBC_GCRS)
        c.representation_type = 'cartesian'

        return np.array([c.x.value, c.y.value, c.z.value])


    def __get_r_sun_unit_GCRS(self,t,loc):
        """
        Returns the cartesian unit vector pointing to the sun from the specified location in the specified instant in
            the Geocentric Celestial Reference Frame (GCRS)
        :param t: current epoch as a Time object <class 'astropy.time.core.Time'>
        :param loc: desired origin location as an Earth Location object <class 'astropy.coordinates.earth.EarthLocation'>
        :return: unit vector pointing to the sun in GCRS [ , , ]
        """
        r_sun_gcrs = get_body('sun', t, loc)
        r_sun_mod = r_sun_gcrs.distance
        r_sun_gcrs.representation_type ='cartesian'
        r_sun_unit = np.array([r_sun_gcrs.x/r_sun_mod, r_sun_gcrs.y/r_sun_mod, r_sun_gcrs.z/r_sun_mod])
        return r_sun_unit

    def __r_SS_SBREF_2_GCRS(self):
        """
        Computes the unit vector in the Geocentric Celestial Reference Frame (GCRS) frame pointing to the sky
            coordinates of the central pixel of the sensor, and stores it in self.__r_SS_SBREF. This method can only run
            if the equivalent vector in the Satellite-Based, Rotating, Ecliptic-Fixed (SBREF) is known and stored in
            self.__r_SS_SBREF.
            NOTE: it does not use the matrices computed by compute_matrices() since the method set_r_SS_SBREF calls it
                before knowing the required attitude
        :return:
        """
        M = self.__get_M_SBREF2GCRS()
        self.__r_SS_GCRS = np.matmul(M, self.__r_SS_SBREF)
        self.__aberration_correction()

    def __aberration_correction(self):
        """
        Applies the necessary correction to the pointing direction due to dilation caused by relativity
        #todo: incomplet; parlar-ho amb Guillem
        :return: N/A
        """

        n = -self.__r_SS_GCRS
        v_sat_geo = self.__orbit.v_geo_gcrs
        earth = get_body_barycentric_posvel('Earth', self.__t)
        vx = earth[1].x.value / (24 * 3600)
        vy = earth[1].y.value / (24 * 3600)
        vz = earth[1].z.value / (24 * 3600)

        v_sat = [v_sat_geo[0]+vx, v_sat_geo[1]+vy, v_sat_geo[2]+vz] #todo: revisar-ho amb les deinicions oficials; no es poden sumar així com així!!

        s = -n + np.cross((1/c)*n,np.cross(v_sat,n)) + (1/c**2)*( np.cross(n*np.dot(n,v_sat),np.cross(v_sat,n) ) + 0.5*np.cross(v_sat,np.cross(n,v_sat)))
        self.__r_SS_GCRS = s

    def __r_SS_GCRS_2_SBREF(self):
        """
        Computes the unit vector in the Satellite-Based, Rotating, Ecliptic-Fixed (SBREF) frame pointing to the sky
            coordinates of the central pixel of the sensor, and stores it in self.__r_SS_SBREF. This method can only run
            if the equivalent vector in the Geocentric Celestial Reference Frame (GCRS) is known and stored in
            self.__r_SS_GCRS.
        :return:
        """
        M = self.__get_M_SBREF2GCRS()
        self.__r_SS_SBREF = np.linalg.inv(M)@self.__r_SS_GCRS

    def __angles_2_r_SS_SBREF(self):
        """
        Calculates the unit vector in the Satellite-Based, Rotating, Ecliptic-Fixed (SBREF) frame pointing to the sky
            point seen in the central pixel of the sensor, given that the attitude and siderostat angles are known and
            stored in self.__yaw, self.__pitch, self.__roll and self.__alpha.
        :return:
        """
        self.__r_SS_SBREF = self.__M_SBBF_SBREF @ self.__M_SBSSF_SBBF @ self.__M_SBFOVF_SBSSF @ [0, 1, 0]

    def __r_SS_SBREF_2_angles(self):
        """
        Defines the necessary attitude (yawm, pitch, roll) and siderostat (alpha) configuration to make the field of
            view axis point to the sky coordinates specified by the unit vector in the Satellite-based, Rotating,
            Ecliptic-Fixed (SBREF) frame stored in the attribute __r_SS_SBREF. The yaw and roll angles are always set
            to zero, and the target attitude is met by the only resulting solution of pitch and alpha. Note that if
            yaw or roll were not set to zero there would be multiple possible solutions, yet out from the nominal
            scanning configuration (which only varies pitch and alpha).
            # todo: afegir l'opció de donar yaw o roll no nuls si l'alpha resultant implica apantallament pel panell solar?

        :return: Prints the resulting alpha and pitch angles [deg], [deg]
        """

        self.__alpha = np.arccos(self.__r_SS_SBREF[1])
        if self.__r_SS_SBREF[2] < 0:
            self.__alpha = 2*pi - self.__alpha
        self.__yaw = 0.0
        self.__pitch = np.arccos(-self.__r_SS_SBREF[0]/np.sin(self.__alpha))
        self.__roll = 0.0
        print(f"alpha = {self.__alpha * 180/pi}")
        print(f"pitch = {self.__pitch * 180/pi}")

    def __check_sun(self):
        """
        Checks if the sun is inside the field of view of the siderostat at the current epoch with the current attitude
            and orbital configuration.
        """
        t = self.__t
        r = self.__orbit.r_geo_inert
        loc = EarthLocation.from_geocentric(r[0], r[1], r[2], u.km)
        r_sun_GCRS = self.__get_r_sun_unit_GCRS(t, loc)
        r_sun_SBFOVF = self.GCRS_2_SBFOVF(r_sun_GCRS)
        if self.__optic.check_if_inframe(r_sun_SBFOVF):
            self.__orbit.print_info()
            raise Exception("The sun is in the field of view!")

    def __check_earth(self):
        """
        Checks if any part of the Earth or its atmosphere (considering an atmospheric height of 100km) is inside the
            field of view of the siderostat at the current epoch with the current attitude and orbital configuration.
        """
        r_E = -self.__orbit.r_geo_inert
        r_E_mod = np.linalg.norm(r_E)
        r_E_u = r_E / r_E_mod
        r_ss = self.__r_SS_GCRS
        h_atm = 100  # [km]
        angle_lim = np.arctan((R+h_atm)/r_E_mod) + max([self.__optic.fov_hor, self.__optic.fov_vert])
        angle = np.abs(np.arccos(np.dot(r_E_u, r_ss)/(np.linalg.norm(r_E_u)*np.linalg.norm(r_ss))))

        if angle < angle_lim:
            #self.__orbit.print_info() # todo: fix bug - contact_gs està buit i no pot imprimir la info quan al primer instant passa això
            #raise Exception("The siderostat is pointing towards the Earth!")
            print("Warning: The siderostat is pointing towards the Earth!")

    def __check_solar_panel(self):
        """
        Checks if the siderostat field of view is pointing to (i.e., blocked by) the solar panel
        """
        if self.__alpha < self.__solar_panel_screening_angle/2 or self.__alpha > (2*pi -self.__solar_panel_screening_angle/2):
            raise Exception("The siderostat is pointing to the solar panel!")

    def __check_attitude_warnings(self):
        """
        Runs all three attitude warning checks (Earth, Sun and Solar panel view blocking)
        :return:
        """
        self.__check_solar_panel()
        self.__check_earth()
        self.__check_sun()

    @property
    def r_SS_SBREF(self):
        """
        :return: The unit vector indicating the direction of the central pixel in the field of view, in the Satellite-
            Based, Rotating, Ecliptic-Fixed frame (SBREF)
        """
        return self.__r_SS_SBREF

    @property
    def angles(self):
        """
        :return: Array with the current three attitude angles, yaw, pitch and roll, in this order [deg, deg, deg]
        """
        return [self.__yaw*180/pi, self.__pitch*180/pi, self.__roll*180/pi, self.__alpha*180/pi]

    @property
    def yaw(self):
        """
        :return: Current Yaw angle [deg]
        """
        return self.__yaw * 180/pi

    @property
    def pitch(self):
        """
        :return: Current Pitch angle [deg]
        """
        return self.__pitch * 180/pi

    @property
    def roll(self):
        """
        :return: Current Roll angle [deg]
        """
        return self.__roll * 180/pi

    @property
    def r_SS_GCRS(self):
        """
        :return: Unit vector in the Geocentric Celestial Reference Frame (GCRS) frame pointing to the sky
            coordinates of the central pixel of the sensor
        """
        return self.__r_SS_GCRS

    @property
    def dec_ICRS(self):
        """
        :return: ICRS declination of the sky point seen in the central pixel of the frame with the current attitude and
            orbital configuration [deg]
        """
        return self.__dec_ICRS

    @property
    def ra_ICRS(self):
        """
        :return: ICRS right ascension of the sky point seen in the central pixel of the frame with the current attitude
            and orbital configuration [deg]
        """
        return self.__ra_ICRS

    #@r_SS_SBREF.setter
    def set_r_SS_SBREF(self, vec):
        """
        Sets the attitude (yaw, pitch, roll) and siderostat pointing (alpha) parameters to make the sky coordinates
            expressed by a unit vector in the Satellite-Based, Rotating, Ecliptic-Fixed (SBREF) frame lay in the center
            of the frame, and warns of Earth, Sun or solar panels screening if vision blocking is found.
        :param vec: unit vector in the SBREF frame with the direction the siderostat is aimed to [ , , ]
        :return: N/A
        """
        self.__r_SS_SBREF = vec
        self.__r_SS_SBREF_2_GCRS()
        self.__r_SS_SBREF_2_angles()
        self.__compute_matrices()
        self.__r_SS_GCRS_2_angles_ICRS()

        self.__check_attitude_warnings()

    #@angles.setter
    def set_angles(self, yaw, pitch, roll, alpha):
        """
        Sets the attitude parametres of the satellite to the given values of yaw, pitch and roll angles, and the
            siderostat position to the given alpha angle. Calculates the sky coordinates to which the siderostat points
            under these conditions (accessible through dec_ICRS() and ra_ICRS() properties) and  warns of Earth, Sun or
            solar panels screening if vision blocking is found.
        :param yaw: desired yaw angle [deg]
        :param pitch: desired pitch angle [deg]
        :param roll: desired roll angle [deg]
        :param alpha: desired siderostat pointing angle [deg]
        :return: N/A
        """
        self.__yaw = yaw * np.pi/180
        self.__pitch = pitch * np.pi/180
        self.__roll = roll * np.pi/180
        self.__alpha = alpha * np.pi/180

        self.__compute_matrices()

        self.__angles_2_r_SS_SBREF()
        self.__r_SS_SBREF_2_GCRS()
        self.__r_SS_GCRS_2_angles_ICRS()

        self.__check_attitude_warnings()

    #@target.setter
    def set_pointing(self, ra_ICRS, dec_ICRS):
        """
        Sets the attitude (yaw, pitch, roll) and siderostat pointing (alpha) parameters to make the sky coordinates
            expressed by a right ascension and a declination in the International Celestial Reference System (ICRS) lay
            in the center of the frame, and warns of Earth, Sun or solar panels screening if vision blocking is found.
        :param ra_ICRS: right ascension of the target point in ICRS coordinates [deg]
        :param dec_ICRS: declination of the target point in ICRS coordinates [deg]
        :return: N/A
        """

        self.__ra_ICRS = ra_ICRS
        self.__dec_ICRS = dec_ICRS

        self.__r_SS_GCRS = self.__angles_ICRS_2_r_GCRS(ra_ICRS, dec_ICRS)
        self.__r_SS_GCRS_2_SBREF()

        self.__r_SS_SBREF_2_angles()
        self.__compute_matrices()

        self.__check_attitude_warnings()

    def __angles_ICRS_2_r_GCRS(self, ra_ICRS, dec_ICRS):
        """
        Converts the sky coordinates given by a right ascension and a declination in the International Celestial
            Reference System (ICRS) to a unit vector in the Geocentric Celestial Reference System (GCRS) pointing to the
            same point.
        :param ra_ICRS: right ascension of a sky point in ICRS coordinates [deg]
        :param dec_ICRS: declination of a sky point in ICRS coordinates [deg]
        :return: r_GCRS: unit vector in GCRS cordinates pointing to the input point
        """

        target_icrs = SkyCoord(ra=ra_ICRS, dec=dec_ICRS, obstime=self.__t, frame='icrs', unit='deg')
        target_gcrs = target_icrs.transform_to('gcrs')
        target_gcrs.representation_type = 'cartesian'
        r_gcrs = np.array([target_gcrs.x.value, target_gcrs.y.value, target_gcrs.z.value])
        target_icrs.representation_type = 'cartesian'

        return r_gcrs


    def __r_SS_GCRS_2_angles_ICRS(self):
        """
        Given the direction where the central point of the field of view is pointing in the form of the unit vector
            r_ss_GCRS in the Geocentric Celestial Reference System (GCRS), computes the sky coordinates of this point in
            the International Celestial Reference System (ICRS). These are stored and accessible through the properties
            dec_ICRS() and ra_ICRS().
        :return: N/A
        """
        coord_GCRS = SkyCoord(x=self.__r_SS_GCRS[0], y=self.__r_SS_GCRS[1], z=self.__r_SS_GCRS[2], obstime=self.__t, frame='gcrs', representation_type='cartesian')
        coord_ICRS = coord_GCRS.transform_to('icrs')
        self.__ra_ICRS = coord_ICRS.ra.value
        self.__dec_ICRS = coord_ICRS.dec.value

    
############################################################################################

    def __compute_matrices(self):
        """
        Computes all the rotation matrices between the different reference systems treated in the calculations, i.e.
            - Satellite-Based, Rotating, Ecliptic-Fixed (SBREF) to Geocentric Celestial (GCRS) and vice-versa;
            - Satellite-Based, Body-Fixed (SBBF) to Satellite-Based, Rotating, Ecliptic-Fixed (SBREF) and vice-versa;
            - Satellite-Based, SideroStat-Fixed (SBSSF) to Satellite-Based, Body-Fixed (SBBF) and vice-versa;
            - Satellite-Based, Field-Of-View-Fixed (SBFOVF) to Satellite-Based, SideroStat-Fixed (SBSSF), and vice-versa
            - Satellite-Based, Field-Of-View-Fixed (SBFOVF) to pqr, and vice-versa
        :return: N/A
        """

        self.__M_SBREF_GCRS = self.__get_M_SBREF2GCRS()
        self.__M_GCRS_SBREF = np.linalg.inv(self.__M_SBREF_GCRS)

        self.__M_SBBF_SBREF = auxf.R1(self.__roll) @ auxf.R2(self.__pitch) @ auxf.R3(self.__yaw)
        self.__M_SBREF_SBBF = np.linalg.inv(self.__M_SBBF_SBREF)

        self.__M_SBSSF_SBBF = auxf.R3(self.__alpha)
        self.__M_SBBF_SBSSF = np.transpose(self.__M_SBSSF_SBBF)

        self.__M_SBFOVF_SBSSF = auxf.R2(self.__alpha)
        self.__M_SBSSF_SBFOVF = np.transpose(self.__M_SBFOVF_SBSSF)

        self.__M_SBFOVF_pqr = np.array([[1,0,0],[0,0,1],[0,-1,0]])
        self.__M_pqr_SBFOVF = np.transpose(self.__M_SBFOVF_pqr)

    def ICRS_2_SBFOVF(self, ra_ICRS, dec_ICRS):
        """
        Converts the sky coordinates given by a right ascension and a declination in the International Celestial
            Reference System (ICRS) to a unit vector in the Satellite-Based, Field-Of-View-Fixed (SBFOVF) frame, based
            on the satellite's current orbital and attitude configuration
        :param ra_ICRS: right ascension of a sky point in ICRS coordinates [deg]
        :param dec_ICRS: declination of a sky point in ICRS coordinates [deg]
        :return: Unit vector in the Satellite-Based, Field-Of-View-Fixed (SBFOVF) frame pointing to the input sky
            coordinates [ , , ]
        """
        
        r_GCRS = self.__angles_ICRS_2_r_GCRS(ra_ICRS, dec_ICRS)
        r_SBFOVF = self.GCRS_2_SBFOVF(r_GCRS)

        return r_SBFOVF

    def ICRS_2_pqr(self, ra_ICRS, dec_ICRS):
        """
        Converts the sky coordinates given by a right ascension and a declination in the International Celestial
            Reference System (ICRS) to a unit vector in the pqr frame, based on the satellite's current orbital and
            attitude configuration
        :param ra_ICRS: right ascension of a sky point in ICRS coordinates [deg]
        :param dec_ICRS: declination of a sky point in ICRS coordinates [deg]
        :return: Unit vector in pqr frame pointing to the input sky coordinates [ , , ]
        """
        r_SBFOVF = self.ICRS_2_SBFOVF(ra_ICRS, dec_ICRS)
        r_pqr = self.SBFOVF_2_pqr(r_SBFOVF)
        return r_pqr

    def SBFOVF_2_pqr(self, r_SBFOVF):
        """
        Converts the sky coordinates given by a unit vector in the Satellite-Based, Field-of-View-Fixed (SBFOVF) frame
            to a unit vector in the pqr frame. Since both reference frames are attached to the optical field of view,
            this transformation is independent of the orbital and attitude parameters.
        :param r_SBFOVF: unit vector in the SBFOVF frame [ , , ]
        :return: r_pqr: unit vector in the pqr frame [ , , ]
        """
        r_pqr = self.__M_SBFOVF_pqr @ r_SBFOVF
        return r_pqr

    def GCRS_2_SBFOVF(self,r_GCRS):
        """
        Converts the sky coordinates given by a cartesian unit vector in the Geocentric Celestial Reference System
            (GCRS) to its coordinates in the Satellite-Based, Field-Of-View-Fixed (SBFOVF) frame
        :param r_GCRS: unit vector in GCRS cartesian coordinates [ , , ]
        :return: r_SBOFOVF unit vector in SBFOVF coordinates [ , , ]
        """
        M = self.__M_SBSSF_SBFOVF @ self.__M_SBBF_SBSSF @ self.__M_SBREF_SBBF @ self.__M_GCRS_SBREF
        r_SBFOVF = M @ r_GCRS

        return r_SBFOVF


    def get_sector(self):
        """
        Gets the folder name with the Gaia star data of the sky sector where the siderostat is currently pointing. At
            the current version it always outputs 'test_sector_2'.
        :return: string with the folder name (sub-folder inside star_catalogs) where the data of the sky sector where
            the siderostat is currently pointing is.
        """
        # todo: definir el sector del cel que toqui segons els arxius del catàleg GAIA
        return 'test_sector_2'


    #  old stuff trash

    # def GCRS_2_SBFOVF(self, r):
    #     M = self.__get_M_SBREF2GCRS()
    #     r_SBREF = np.matmul(np.linalg.inv(M), r)
    #
    #     M = np.transpose(auxf.R1(self.__roll) @ auxf.R1(self.__pitch) @ auxf.R1(self.__yaw))
    #     r_SBBF =  np.matmul(M, r_SBREF)
    #
    #     M = auxf.R1(self.__alpha)
    #     r_SBSSF = np.matmul(M, r_SBBF)
    #
    #     M = auxf.R1(self.__alpha)
    #     r_SBFOVF = np.matmul(M, r_SBSSF)
    #
    #     return r_SBFOVF


    # def __get_k_ICRS_GCRS(self):
    #     k_icrs = SkyCoord(ra=0, dec=90, obstime=self.__t, frame='icrs', unit='deg')
    #     k_icrs_gcrs = k_icrs.transform_to('gcrs')
    #     k_icrs_gcrs.representation_type = 'cartesian'
    #     return [k_icrs_gcrs.x, k_icrs_gcrs.y, k_icrs_gcrs.z]

    # def __get_r_sun_GCRS(self):
    #     loc = EarthLocation.from_geocentric(self.__orbit.__r_geo_rot[0], self.__orbit.__r_geo_rot[1], self.__orbit.__r_geo_rot[2], u.km)
    #     r_sun_gcrs = get_body('sun', self.__t, loc)
    #     r_sun_mod = r_sun_gcrs.distance
    #     r_sun_gcrs.representation_type = 'cartesian'
    #     return [(r_sun_gcrs.x / r_sun_mod), (r_sun_gcrs.y / r_sun_mod), (r_sun_gcrs.z / r_sun_mod)]