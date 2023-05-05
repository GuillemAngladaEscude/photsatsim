import numpy as np

from satellite_class import Satellite
import aux_functions as auxf
from astropy.time import Time
from astropy import units as u

import astropy.coordinates
from astropy.coordinates import solar_system_ephemeris
solar_system_ephemeris.set('jpl') #'de432s'


from Image import Image
from optic_class import Optic
from orbit_class import Orbit


#############################################################################
# Definim els paràmetres orbitals de PhotSat
zp = 500  # altura del perigeu [km]
za = 550  # altura de l'apogeu [km]
i = auxf.find_i_SSO(zp,za)  # inclinació necessària per garantir òrbita heliosíncrona [deg]

# Definim el temps inicial de la simulació - ara posat al solstici d'estiu
t0 = Time('2023-06-24 12:00:00.000', scale='tcb')

# objecte òrbita
orbit_photsat = Orbit.from_elements(t0, zp, za, i, raan=0, theta_0=90, aop=0)
# afegim les estacions del Montsec i, opcionalment, Svalbard
# això serveix afegir l'opció de no generar imatges quan hi hagi contacte amb una estació de terra per simular que s'estan enviant dades
orbit_photsat.add_ground_station("Montsec")  # a cada instant de la propagació hi haurà un 1 o un 0 a photsat.orbit.gs_contact[0]
#orbit_photsat.add_ground_station("Svalbard")  # a cada instant de la propagació hi haurà un 1 o un 0 a photsat.orbit.gs_contact[1]


# variables òptiques a definir
N_pixels_hor = 150  # nombre de píxels per costat del sensor
N_pixels_vert = 150
fov = 6  # field of view [deg]
focal_length = 10  # focal length [cm]

max_mag = 9  # magnitud màxima a què limitem el catàleg a carregar

# variables òptiques dependents:
d_hor = 2*focal_length*np.tan(fov/2*np.pi/180)   # sensor horizontal size [cm]
d_vert = 2*focal_length*np.tan(fov/2*np.pi/180)  # sensor vertical size [cm]

# objecte òptica
optic_photsat = Optic(N_pixels_hor, N_pixels_vert, d_hor, d_vert, focal_length, max_mag)

# objecte satèl·lit
photsat = Satellite(optic_photsat, orbit_photsat)


solar_panel_screening_angle = 90
s_AL = 6  # salt along_track entre dues imatges consecutives [deg]
s_C = 5   # salt cross_track entre dues imatges consecutives [deg]
# nota: along-track i cross-track estan definits respecte el cercle màxim resseguit a l'esfera celeste i no respecte
# l'òrbita de PhotSat

pics_orbit = (360 - solar_panel_screening_angle - fov)/s_AL + 1  # imatges generades per cada òrbita
orbits = 180/s_C  # òrbites necessàries per cobrir tota l'esfera celeste
# todo: afegir warnings de nombres enters per pics_orbit i orbits


delta_t = (photsat.orbit.T/60)/pics_orbit  # temps d'integració resultant per a cada imatge
print(f"delta_t = {delta_t}")

i = 0  # comptador d'imatges generades
for j in range(round(orbits)):
    # actualitzem el pitch, que es manté constant durant tota l'òrbita
    pitch = j*s_C
    print(f"pitch = {pitch}")
    for k in range(round(pics_orbit)):
        print(f"Picture number: {i}")
        # actualitzem l'angle alpha, que va augmentant progressivament entre foto i foto
        alpha = (solar_panel_screening_angle + fov)/2 + k*s_AL
        print(f"alpha = {alpha}")

        # propaguem la posició de PhotSat a l'instant actual
        photsat.orbit.propagate(delta_t*i)
        photsat.att.set_angles(0, pitch, 0, alpha)
        print(f"Pointing to ra = {photsat.att.ra_ICRS}, dec = {photsat.att.dec_ICRS}")

        stars_in_frame = photsat.get_stars_in_frame()

        # Generem la imatge
        image = Image(N_pixels_vert, N_pixels_hor)
        image.placeStar(stars_in_frame)

        # Desem la imatge (nom actual: FullScan_Image<i>)
        image.SaveToFile(i)

        i = i+1
