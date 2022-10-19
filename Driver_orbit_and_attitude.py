import numpy as np

from satellite_class import Satellite
import aux_functions as auxf
from astropy.time import Time
from astropy import units as u

import astropy.coordinates
from astropy.coordinates import SkyCoord
from astropy.coordinates import solar_system_ephemeris, EarthLocation, CartesianRepresentation
from astropy.coordinates import get_body, get_sun, get_body_barycentric, get_body_barycentric_posvel
solar_system_ephemeris.set('jpl') #'de432s'


from Image import Image
from optic_class import Optic
from catalog_class import Catalog
from orbit_class import Orbit




#file = open("/Users/llucr/PycharmProjects/PhotSat/data_stars_upto_mag_5.csv")
#type(file)

#############################################################################
# obrim csv del catàleg GAIA
max_mag = 14  # magnitud màxima a què limitem el catàleg a carregar
# sector = 'test_sector_2'
# stardata = Catalog.load(max_mag,sector)
# print("Star catalog loaded successfully!")
#############################################################################

#############################################################################
# Definim els paràmetres orbitals de PhotSat
zp = 550 # altura del perigeu [km]
za = 1000 # altura de l'apogeu [km]
i = auxf.find_i_SSO(zp,za) #inclinació necessària per garantir òrbita heliosíncrona [deg]
t0 = Time('2023-06-24 12:00:00.000', scale='tcb')  # temps inicial de la simulació # tcb/tai/ut1?

orbit_photsat = Orbit.from_elements(t0, zp, za, i, raan = 0, theta_0 = 20, aop = 0)
orbit_photsat.add_ground_station("Montsec")
orbit_photsat.add_ground_station("Svalbard")

# r = np.array([6024.52465821, -512.9644707, 3440.22769172])
# v = np.array([-3.73427845, -0.98852292, 6.62958961])
# orbit_photsat = Orbit.from_rv(t0,r,v)

N_pixels_hor = 100  # nombre de píxels per costat del sensor (quadrat) []
N_pixels_vert = 100
fov = 6  # field of view [deg]
focal_length = 10  # focal length [cm]
d_hor = 2*focal_length*np.tan(fov/2*np.pi/180)   # sensor horizontal size [cm]
d_vert = 2*focal_length*np.tan(fov/2*np.pi/180)  # sensor vertical size [cm]
optic_photsat = Optic(N_pixels_hor, N_pixels_vert, d_hor, d_vert, focal_length, max_mag)
#optic_photsat = Optic(N_pixels, fov, focal_length, max_mag)

photsat = Satellite(optic_photsat, orbit_photsat)

#############################################################################

#############################################################################
# definim les condicions a l'instant que volem avaluar
delta_t = 125 #temps transcorregut entre t0 i l'instant en què ens trobem [min]
t = t0 + astropy.time.core.TimeDelta(delta_t*u.min, scale='tcb') # ut1 ?
photsat.orbit.propagate(delta_t)
photsat.orbit.print_info()

#############################################################################

#############################################################################
# definim el target: volem que el Photsat apunti a l'estel i del catàleg
# i = 1 #índex del catàleg reduït que hem importat
# ra_ICRS = stardata[i][2]
# dec_ICRS = stardata[i][3]

# ra_ICRS = 47
# dec_ICRS = 7
ra_ICRS = 81
dec_ICRS = -69
print(f"Configurem Photsat perquè apunti a les coordenades ICRS ra = {ra_ICRS} deg; dec = {dec_ICRS} deg:")
# calculem angles alpha i pitch necessaris perquè les coordenades target quedi al centre del frame
photsat.att.set_pointing(ra_ICRS, dec_ICRS)
#############################################################################

#############################################################################
# construïm el catàleg de tots els estels que quedin dins del frame en la posició actual
stars_in_frame = photsat.get_stars_in_frame()
#############################################################################

#############################################################################

# Creem una nova imatge
image = Image(N_pixels_vert, N_pixels_hor)

#afegim els estels de la llista stars_in_frame a la imatge:
image.placeStar(stars_in_frame)

image.plotImage()
print("Process finished")