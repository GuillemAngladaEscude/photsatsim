import csv

class Catalog:

    def __init__(self):
        pass

    @classmethod
    def load(self, max_mag, sector):
        # obrim csv del catàleg GAIA
        import os
        THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
        my_file = os.path.join(THIS_FOLDER, 'star_catalogs', sector,
                               'data_stars_upto_mag_' + str(max_mag) + '.csv')  # ara mateix obrim la versió filtrada a mag<=5
        file = open(my_file)
        csvreader = csv.reader(file)
        #header = []
        header = next(csvreader)
        stardata = []
        for row in csvreader:
            stardata.append(row)
        file.close()

        return stardata