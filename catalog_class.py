import csv
from attitude_class import Attitude
from optic_class import Optic
from astroquery.gaia import Gaia
import numpy as np

class Catalog:
    """
    A class used to comprise the methods to load the catalog of stars visible by PhotSat given its optical instrument
    and its current attitude.

    This class uses the modules numpy, astroquery.gaia and csv.

    Attributes
    ----------


    Methods
    -------
    load(optic, attitude, method):
        Given the optic and attitude objects for PhotSat at the current instant and a method of catalog loading, it
        outputs the list of stars visible in the frame. The resulting array is n rows x 4 columns, each of which
        containing the star ID, ICRS right ascension [deg], ICRS declination [deg] and magnitude [mean phot. g]

    """

    def __init__(self):
        pass

    @classmethod
    def load(cls, optic: Optic, att: Attitude, method):
        """"
        Computes the list of stars inside the frame at the current instant. If the method is set to
        'direct_gaia_download', a request to the GAIA catalog is sent online with all the stars at a radial distance of
        max([optic.fov_hor, optic.fov_vert])/sqrt(2) of the central point of the frame. Otherwise, if the method is set
        to 'sector_load', a csv file with the corresponding stars is loaded. This assumes the whole GAIA catalog has
        been divided in different csv files locally available, each file sorresponding to a given sector of the sky and
        a maximum star magnitude. This division is currently not implemented and thus the method get_sector() inside the
        attitude class always returns 'test_sector_2'.

        :param optic: Optic object corresponding to the Optic system on board PhotSat
        :param att: Attitude object corresponding to the attitude of the spacecraft at the current instant

        :return: stardata: double array with n rows x 4 columns, each of which containing the star ID, ICRS right
            ascension [deg], ICRS declination [deg] and magnitude [mean phot. g]
        """

        max_mag = optic.mag_max

        if method == 'direct_gaia_download':

            fov = max([optic.fov_hor, optic.fov_vert])*180/np.pi
            ra = att.ra_ICRS
            dec = att.dec_ICRS

            job = Gaia.launch_job_async(f"""
            SELECT gaia_source.source_id, gaia_source.ra,gaia_source.dec,gaia_source.phot_g_mean_mag,
            DISTANCE(
               POINT({ra}, {dec}),
               POINT(ra, dec)) AS ang_sep
            FROM gaiadr3.gaia_source
            WHERE 1 = CONTAINS(
               POINT({ra}, {dec}),
               CIRCLE(ra, dec, {fov / np.sqrt(2)}))
            AND phot_g_mean_mag < {max_mag}
            ORDER BY ang_sep ASC
            """, dump_to_file=False)

            r = job.get_results()
            id_vec = np.array(r['source_id'].data)
            ra_vec = np.array(r['ra'].data)
            dec_vec = np.array(r['dec'].data)
            mag_vec = np.array(r['phot_g_mean_mag'].data)

            stardata = []
            for i in range(len(r)):
                stardata.append([id_vec[i], ra_vec[i], dec_vec[i], mag_vec[i]])



        elif method == 'sector_load':
            # obrim csv del catàleg GAIA
            sector = att.get_sector()
            import os
            THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
            my_file = os.path.join(THIS_FOLDER, 'star_catalogs', sector,
                               'data_stars_upto_mag_' + str(max_mag) + '.csv')
            file = open(my_file)
            csvreader = csv.reader(file)
            #header = []
            header = next(csvreader)
            stardata = []
            for row in csvreader:
                stardata.append(row)
            file.close()

        else:
            raise Exception(f"Method {method} is not recognized as a valid method for star catalog obtention!")

        return stardata