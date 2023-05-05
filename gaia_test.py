import astropy.units as u
import matplotlib
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astroquery.gaia import Gaia

#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np

# Suppress warnings. Comment this out if you wish to see the warning messages
# import warnings
# warnings.filterwarnings('ignore')

from astroquery.gaia import Gaia
# tables = Gaia.load_tables(only_names=True)
# for table in (tables):
#     print (table.get_qualified_name())

# job = Gaia.launch_job_async("SELECT * \
#     FROM gaiadr1.gaia_source \
#     WHERE CONTAINS(POINT(gaiadr1.gaia_source.ra,gaiadr1.gaia_source.dec),CIRCLE(56.75,24.1167,2))=1;" \
#                                 , dump_to_file=False)

mag_max = 12
ra = 0
dec = 90
fov = 6

job = Gaia.launch_job_async(f"""
SELECT gaia_source.ra,gaia_source.dec,gaia_source.phot_g_mean_mag,
DISTANCE(
   POINT({ra}, {dec}),
   POINT(ra, dec)) AS ang_sep
FROM gaiadr3.gaia_source
WHERE 1 = CONTAINS(
   POINT({ra}, {dec}),
   CIRCLE(ra, dec, {fov/np.sqrt(2)}))
AND phot_g_mean_mag < {mag_max}
ORDER BY ang_sep ASC
""", dump_to_file=False)

r = job.get_results()
ra_vec = (np.array(r['ra'].data))
print(ra_vec)
print(r)