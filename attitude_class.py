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

import aux_functions as auxf

c = 299792.458  # km/s

class Attitude:

    def __init__(self, orbit: Orbit):

        #parking values
        self.__yaw = 0
        self.__pitch = 0
        self.__roll = 0
        self.__alpha = 180

        self.__orbit = orbit
        self.__t = orbit.t

        #self.set_angles(0,0,0,180)
        self.__angles_2_r_SS_SBREF()
        self.__r_SS_SBREF_2_GCRS()
        self.__r_SS_GCRS_2_angles_ICRS()

    def __get_M_SBREF2GCRS(self):

        t = self.__orbit.t_last_perigee
        M = np.zeros((3,3))

        EM_barycenter_GCRS = get_body('earth-moon-barycenter', t)  # GCRS: Geocentric Celestial Reference Frame
        EM_barycenter_GCRS.representation_type = 'cartesian'
        r = [EM_barycenter_GCRS.x.value, EM_barycenter_GCRS.y.value, EM_barycenter_GCRS.z.value]
        loc = EarthLocation.from_geocentric(r[0], r[1], r[2], u.km)
        a = EarthLocation.get_gcrs_posvel(loc, t)
        r_EMBC_GCRS = a[0]  # EMBC: Earth-Moon Barycenter
        v_EMBC_GCRS = a[1]

        BME_frame = astropy.coordinates.builtin_frames.BarycentricMeanEcliptic()

        a = SkyCoord(lon=0, lat=90, obstime=t, frame=BME_frame, unit='deg')
        b = a.transform_to('gcrs')

        c = astropy.coordinates.GCRS(ra=b.ra, dec=b.dec, obstime=t, obsgeoloc=r_EMBC_GCRS, obsgeovel=v_EMBC_GCRS)
        c.representation_type = 'cartesian'

        k_ecl_GCRS = np.array([c.x.value, c.y.value, c.z.value])

        M[0,2] = k_ecl_GCRS[0]
        M[1,2] = k_ecl_GCRS[1]
        M[2,2] = k_ecl_GCRS[2]

        r_sun_gcrs = get_body('sun', t, loc)
        r_sun_mod = r_sun_gcrs.distance
        r_sun_gcrs.representation_type ='cartesian'
        r_sun_unit = [r_sun_gcrs.x/r_sun_mod, r_sun_gcrs.y/r_sun_mod, r_sun_gcrs.z/r_sun_mod]
        M[0,1] = r_sun_unit[0]
        M[1,1] = r_sun_unit[1]
        M[2,1] = r_sun_unit[2]

        x = np.cross(r_sun_unit, k_ecl_GCRS)
        M[0,0] = x[0]
        M[1,0] = x[1]
        M[2,0] = x[2]

        #r_sun_gcrs = get_body('sun', self._Orbit__t, loc)
        #print(f"SBREF base determinant: {np.linalg.det(M)}")
        #print(f"M_SBREF_GCRS: {M}")
        #print(f"M_GCRS_SBREF: {np.linalg.inv(M)}")
        return M

    def __r_SS_SBREF_2_GCRS(self): #TODO: canviar a target i afegir transformacions  a ICRS?
        M = self.__get_M_SBREF2GCRS()
        self.__r_SS_GCRS = np.matmul(M,self.__r_SS_SBREF)
        self.__aberration_correction()

    def __aberration_correction(self):

        n = -self.__r_SS_GCRS
        v_sat_geo = self.__orbit.v_geo_gcrs
        earth = get_body_barycentric_posvel('Earth', self.__t)
        vx = earth[1].x.value / (24 * 3600)
        vy = earth[1].y.value / (24 * 3600)
        vz = earth[1].z.value / (24 * 3600)

        v_sat = [v_sat_geo[0]+vx, v_sat_geo[1]+vy, v_sat_geo[2]+vz] #todo: revisar-ho amb les deinicions ofiials; no es poden sumar així com així!!

        s = -n + np.cross((1/c)*n,np.cross(v_sat,n)) + (1/c**2)*( np.cross(n*np.dot(n,v_sat),np.cross(v_sat,n) ) + 0.5*np.cross(v_sat,np.cross(n,v_sat)))
        self.__r_SS_GCRS = s

    def __r_SS_GCRS_2_SBREF(self): #todo: revisar
        M = self.__get_M_SBREF2GCRS()
        #self.__r_SS_GCRS = np.matmul(np.linalg.inv(M), self.__r_SS_GCRS)
        self.__r_SS_SBREF = np.linalg.inv(M)@self.__r_SS_GCRS

    def __angles_2_r_SS_SBREF(self):
        angle_SS = self.__alpha*pi/180
        r_ss_0 = np.matmul(auxf.R3(angle_SS),[0, 1, 0])
        yaw = self.__yaw*pi/180
        pitch = self.__pitch*pi/180
        roll = self.__roll*pi/180
        M = np.transpose(auxf.R1(roll)@auxf.R2(pitch)@auxf.R3(yaw))
        self.__r_SS_SBREF = M @ r_ss_0

    def __r_SS_SBREF_2_angles(self):

        self.__alpha = np.arccos(self.__r_SS_SBREF[1])
        if self.__r_SS_SBREF[2]<0:
            self.__alpha = 2*pi - self.__alpha
        self.__yaw = 0.0
        self.__pitch = np.arccos(-self.__r_SS_SBREF[0]/np.sin(self.__alpha))
        self.__roll = 0.0
        #todo: posar warnings d'apantallament per sol, terra, panell solar
        print(f"alpha = {self.__alpha * 180/pi}")
        print(f"pitch = {self.__pitch * 180/pi}")

    @property
    def r_SS_SBREF(self):
        return self.__r_SS_SBREF

    @property
    def angles(self):
        return [self.__yaw*180/pi, self.__pitch*180/pi, self.__roll*180/pi, self.__alpha*180/pi]

    @property
    def yaw(self):
        return self.__yaw * 180/pi

    @property
    def pitch(self):
        return self.__pitch * 180/pi

    @property
    def roll(self):
        return self.__roll * 180/pi

    @property
    def r_SS_GCRS(self):
        return self.__r_SS_GCRS

    @property
    def dec_ICRS(self):
        return self.__dec_ICRS

    @property
    def ra_ICRS(self):
        return self.__ra_ICRS

    #@r_SS_SBREF.setter
    def set_r_SS_SBREF(self, vec):
        self.__r_SS_SBREF = vec
        self.__r_SS_SBREF_2_GCRS()
        self.__r_SS_SBREF_2_angles()
        self.__compute_matrices()
        self.__r_SS_GCRS_2_angles_ICRS()

    #@angles.setter
    def set_angles(self,yaw,pitch,roll,alpha):
        self.__yaw = yaw
        self.__pitch = pitch
        self.__roll = roll
        self.__alpha = alpha

        self.__compute_matrices()

        self.__angles_2_r_SS_SBREF()
        self.__r_SS_SBREF_2_GCRS()
        self.__r_SS_GCRS_2_angles_ICRS()

    #@target.setter
    def set_pointing(self, ra_ICRS, dec_ICRS):

        self.__ra_ICRS = ra_ICRS
        self.__dec_ICRS = dec_ICRS

        self.__r_SS_GCRS = self.__angles_ICRS_2_r_GCRS(ra_ICRS, dec_ICRS)
        #print(f"r_GCRS: {self.__r_SS_GCRS}")
        self.__r_SS_GCRS_2_SBREF()

        #print(f"r_SBREF: {self.__r_SS_SBREF}")
        self.__r_SS_SBREF_2_angles()
        self.__compute_matrices()

    def __angles_ICRS_2_r_GCRS(self, ra_ICRS, dec_ICRS):
        target_icrs = SkyCoord(ra=ra_ICRS, dec=dec_ICRS, obstime=self.__t, frame='icrs', unit='deg')
        target_gcrs = target_icrs.transform_to('gcrs')
        target_gcrs.representation_type = 'cartesian'
        r_gcrs = np.array([target_gcrs.x.value, target_gcrs.y.value, target_gcrs.z.value])

        target_icrs.representation_type = 'cartesian'
        #print(f"r_ICRS: {target_icrs}")
        return r_gcrs


    def __r_SS_GCRS_2_angles_ICRS(self):
        coord_GCRS = SkyCoord(x=self.__r_SS_GCRS[0], y=self.__r_SS_GCRS[1], z=self.__r_SS_GCRS[2], obstime=self.__t, frame='gcrs', representation_type='cartesian')
        coord_ICRS = coord_GCRS.transform_to('icrs')
        self.__ra_ICRS = coord_ICRS.ra.value
        self.__dec_ICRS = coord_ICRS.dec.value

    
############################################################################################

    def __compute_matrices(self):

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
        
        r_GCRS = self.__angles_ICRS_2_r_GCRS(ra_ICRS, dec_ICRS)
        M = self.__M_SBSSF_SBFOVF @ self.__M_SBBF_SBSSF @ self.__M_SBREF_SBBF @ self.__M_GCRS_SBREF
        r_SBFOVF = M @ r_GCRS

        return r_SBFOVF

    def ICRS_2_pqr(self, ra_ICRS, dec_ICRS):
        r_SBFOVF = self.ICRS_2_SBFOVF(ra_ICRS, dec_ICRS)
        r_pqr = self.__M_SBFOVF_pqr @ r_SBFOVF
        return r_pqr

    def SBFOVF_2_pqr(self, r_SBFOVF):
        r_pqr = self.__M_SBFOVF_pqr @ r_SBFOVF
        return r_pqr


    def get_sector(self):
        # todo: definir el sector del cel que toqui segons els arxius del catàleg GAIA
        return 'test_sector_2'





    #  old stuff trash

    # def GCRS_2_SBFOVF(self, r):
    #     M = self.__get_M_SBREF2GCRS()
    #     r_SBREF = np.matmul(np.linalg.inv(M), r)
    #
    #     M = np.transpose(auxf.R1(self.__roll) @ auxf.R1(self.__pitch) @ auxf.R1(self.__yaw)) #todo: revisar; agrupar en un mètode per no duplicar codi
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