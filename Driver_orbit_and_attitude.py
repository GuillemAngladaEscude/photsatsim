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


#file = open("/Users/llucr/PycharmProjects/PhotSat/data_stars_upto_mag_5.csv")
#type(file)

#############################################################################
# obrim csv del catàleg GAIA
max_mag = 8  # magnitud màxima a què limitem el catàleg a carregar
sector = 'test_sector'
stardata = Catalog.load(max_mag,sector)
print("Star catalog loaded successfully!")
#############################################################################

#############################################################################
# Definim els paràmetres orbitals de PhotSat
zp = 550 # altura del perigeu [km]
za = 1000 # altura de l'apogeu [km]
i = auxf.find_i_SSO(zp,za) #inclinació necessària per garantir òrbita heliosíncrona [deg]
t0 = Time('2023-06-30 00:00:00.000', scale='tcb') # temps inicial de la simulació # tcb o tai?

N_pixels = 100 # nombre de píxels per costat del sensor (quadrat)
fov = 6 # field of view [deg]
optic_photsat = Optic(N_pixels, fov)
photsat = Satellite(optic_photsat,t0,zp,za,i,0,0,0)
#############################################################################

#############################################################################
# definim les condicions a l'instant que volem avaluar
delta_t = 140 #temps transcorregut entre t0 i l'instant en què ens trobem [min]
t = t0 + astropy.time.core.TimeDelta(delta_t*u.min, scale='tcb') # ut1 ?
photsat.orbit.propagate(delta_t)
photsat.orbit.print_elements()
#############################################################################

#############################################################################
# definim el target: volem que el Photsat apunti a l'estel i del catàleg
i = 50 #índex del catàleg reduït que hem importat
ra_ICRS = stardata[i][2]
dec_ICRS = stardata[i][3]
print(f"Configurem Photsat perquè apunti a les coordenades ICRS ra = {ra_ICRS} deg; dec = {dec_ICRS} deg:")
photsat.att.set_pointing(dec_ICRS, ra_ICRS) #calculem angles alpha i pitch necessaris perquè l'estel target quedi al centre del frame
#############################################################################

#############################################################################
# construïm el catàleg de tots els estels que quedin dins del frame en la posició actual
stars_in_frame = photsat.get_stars_in_frame(stardata)
print("stars in frame successfully obtained")
#############################################################################

#############################################################################

# Creem una nova imatge
image = Image(N_pixels,N_pixels)

#afegim els estels de la llista stars_in_frame a la imatge:
image.placeStar(stars_in_frame)

image.plotImage()
print("Process finished")