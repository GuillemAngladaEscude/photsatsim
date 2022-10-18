from orbit_class import Orbit
from attitude_class import Attitude
from catalog_class import Catalog
from optic_class import Optic

from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body, get_sun, get_body_barycentric, get_body_barycentric_posvel
solar_system_ephemeris.set('jpl') #'de432s'


from tqdm import tqdm

class Satellite:


    def __init__(self, optic: Optic, orbit: Orbit):

        #self.orbit = Orbit(t0, zp, za, i, raan, theta_0, aop)
        self.orbit = orbit
        self.optic = optic
        self.att = Attitude(self.orbit)


    def get_stars_in_frame(self): # aquest mètode va aquí perquè depèn tant d'attitude com d'òptica

        print("Obtaining stars in frame...")
        sector = self.att.get_sector()
        stardata = Catalog.load(self.optic.mag_max, sector)

        stars_in_frame = []
        for star in tqdm(stardata, colour="WHITE"):

            # coordenades ICRS de cada estel:
            ra_ICRS = star[2]
            dec_ICRS = star[3]

            # vector que apunta a l'estel en la base Satellite Based Field Of View Fixed:
            r_SBFOVF = self.att.ICRS_2_SBFOVF(ra_ICRS, dec_ICRS)

            # coordenades de l'estel en píxels dins del frame (new):
            if self.optic.check_if_inframe(r_SBFOVF):
                r_pqr = self.att.SBFOVF_2_pqr(r_SBFOVF)
                u, v = self.optic.project_pqr_2_uv(r_pqr)
                x, y = self.optic.uv_2_xy(u, v)

                #frame_coord = self.optic.project_SBFOVF_2_sensor(r_SBFOVF, 'linear') #aquest era per quan no hi havia
                # el mètode xy_2_xy_TL a la classe Image
                frame_coord = [x, y]
                new_row = [frame_coord[0], frame_coord[1],
                           float(star[5])]  # star[5] dona la magnitud aparent; canviar-ho a star[4] si es vol el flux
                stars_in_frame.append(new_row)

        print("Stars in frame successfully obtained!")
        return stars_in_frame

