import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
# import matplotlib.scale.LogScale as logscale
# comentary

class Image:

    def __init__(self, npx, npy, stampSize = 5, libname='none'):
        """This class must be initialized introducing 'npx' & 'npy' which are the number of pixels of the image
        on each side. Also, we can introduce 'stampSize' or 'libname'. If 'stampSize' is introduced, a library
        of stamps with that number of pixels on each side is selected (the default is stamps of 5 pixels on each side).
        If 'libname' is introduced, the image will be generated using the library of stamps that has that name
        (this needs only to be used if we want to work with a library that is not generated with 'Driver_library.py')."""

        # number of pixels of the image
        self.__npx = npx
        self.__npy = npy

        if libname != 'none':
            # if a library name ('libname') is introduced, we call a function that loads a library with this name
            self.__initLibraryWithLibraryName(libname)
        else:
            # if not, we load a library of stamps with the same number of pixels as 'stampSize'
            self.__initLibraryWithStampSize(stampSize)

        # generation of an image with npx * npy pixels, each of them initially with zero flux
        self.__array = np.zeros((self.__npx, self.__npy))

    # set of the number of pixels of the stamp
    def __setPixelsStamp(self, npx_stamp_input, npy_stamp_input):
        """This function sets the number of pixels that the stamps
         have on each side"""
        self.__npx_stamp = npx_stamp_input
        self.__npy_stamp = npy_stamp_input
        # The stamp must have an odd number of pixels at each side. If this does not happen, an error occurs
        if self.__npx_stamp % 2 == 0 or self.__npy_stamp % 2 == 0:
            raise ValueError('The stamp must have an odd number of pixels at each side.')

    def __loadLibrary(self, libraryname):
        """This function loads a file corresponding to a library of stamps
        with the name 'libraryname'. This library is used to generate the image."""
        self.__stars = np.load(libraryname)

    def __initLibraryWithStampSize(self, stampSize):
        """This function initializes a library with the stamps of the size of
        'stampSize', which shall be an odd integer (1,3,5,7,9...)."""
        libname = 'Library_'+str(stampSize)+"x"+str(stampSize)+'.npy'
        self.__setPixelsStamp(stampSize, stampSize)
        self.__loadLibrary(libname)

    def __initLibraryWithLibraryName(self, libname_input):
        """This function initializes a library with the name introduced in
        'libname_input'. Then obtains the number of pixels of the stamps of the library
        (supposing that have the same number of pixels on each side). Finally, sets the number of
        pixels of the stamps."""
        self.__loadLibrary(libname_input)
        stampSize = self.__stars[0][0]
        self.__setPixelsStamp(stampSize, stampSize)

    def __mag_to_flux(self, magnitude):
        """This function receives an apparent magnitude and returns the corresponding
        radiant flux (normalized to be 1 for magnitude = 10)"""

        L_sun = 3.828 * 10 ** 33  # solar luminosity: [erg * s^-1] == radiant flux * surface
        L_vega = 37 * L_sun  # luminosity of Vega

        m_vega = 0.03  # apparent magnitude of Vega
        m_norm = 10  # normalization magnitude: the flux corresponding to this magnitude will be normalized to 1

        F_norm = L_vega / 100 ** (
                    (float(m_norm) - m_vega) / 5.0)  # flux corresponding to m_norm that will be normalized to 1
        F_vega_norm = L_vega / F_norm  # flux of Vega after the normalization

        flux = F_vega_norm / 100 ** (
                    (float(magnitude) - 0.03) / 5.0)  # normalized flux corresponding to the introduced 'magnitude'

        return flux

    def get_array(self):
        """This function returns the saved in memory image (modifying this array means modify the saved
        in memory array = image)."""
        return self.__array

    def get_arrayCopy(self):
        """This function generates a copy of the image (which is a class variable)."""
        copy = self.__array.copy()
        return copy

    def test(self, number_partitions):
        """This function tests that the flux inside the image remains constant
        while we place a star in different points of a pixel. 'number_partitions' corresponds
        to the number of partitions that we want to do of the pixel for the test."""

        pixel_x, pixel_y = 4.357, 5

        dx, dy = 1.0 / number_partitions, 1.0 / number_partitions
        # plots = [[0.] * number_partitions for i in range(number_partitions)]   # np.zeros(number_partitions * number_partitions)
        results = np.zeros(number_partitions)  # , number_partitions))
        positions = np.zeros(number_partitions)

        # for i in range(number_partitions):
        for j in range(number_partitions):  # todo: pensar si volem nomes una linia de pixel o tot. quin resultat volem exactament?

            flux_sum = 0.0

            # position_x = pixel_x + i * dx
            position_y = pixel_y + j * dy
            test_catalog = [[pixel_x, position_y, 1000.0]]

            positions[j] = position_y

            test_image = Image(self.__npx, self.__npy)
            test_image.placeStar(test_catalog)
            # plots[i][j] = test_image.__array

            for I in range(self.__npx):
                for J in range(self.__npy):
                    flux_sum += test_image.__array[I][J]

            results[j] = flux_sum

        final_result = np.matrix.flatten(results)

        maximum = max(final_result)
        minimum = min(final_result)

        main_difference = maximum - minimum
        print(main_difference, maximum, minimum)

        # print(positions)

        plt.plot(positions, final_result)
        plt.xlabel('pixel position')
        plt.ylabel('flux')
        plt.show()

    def placeStar(self, catalog):
        """This function places each star of 'catalog' which must have the shape
        catalog = [[position pixel for x of star 1, position pixel for y of star 1, apparent magnitude of star 1], [idem for star 2], [idem for star 3], ...]"""

        for k in range(len(catalog)):

            # we assign each value of the catalog to its corresponding variable
            x, y, flux = catalog[k][0], catalog[k][1], self.__mag_to_flux(catalog[k][2])

            # we calculate the pixel and the partition of the pixel in which is the star
            pix_x1, par_x = int(x // 1.0), x % 1.0
            pix_y1, par_y = int(y // 1.0), y % 1.0

            # we calculate the number of partitions that have the pixels of the stamps (the precision of the library)
            partitions_x_stamp = len(self.__stars)
            partitions_y_stamp = len(self.__stars[0])

            # with the number of partitions we calculate the length of each partition
            dx_partitions = 1.0 / float(partitions_x_stamp)
            dy_partitions = 1.0 / float(partitions_y_stamp)

            # with this loop we identify between which stamps is our star (for the bilinear interpolation)
            for i in range(partitions_x_stamp + 1):
                if par_x <= (i * dx_partitions):
                    break

            # we assign the result of the loop to these variables in order to use it later to obtain the correct stamps
            par_x1 = i - 1
            par_x2 = i

            # we calculate these variables because are necessary to compute the bilinear interpolation
            real_partition_x1 = par_x1 * dx_partitions
            real_partition_x2 = par_x2 * dx_partitions

            # now we make sure that if the star is in the 0.0 spot, the stamp is correctly assigned
            if par_x1 == -1:
                par_x1 = 0
                par_x2 = 1

            # idem for y (loop to obtain the correct stamps)
            for j in range(partitions_y_stamp + 1):
                if par_y <= (j * dy_partitions):
                    break

            par_y1 = j - 1
            par_y2 = j

            real_partition_y1 = par_y1 * dy_partitions
            real_partition_y2 = par_y2 * dy_partitions

            if par_y1 == -1:
                par_y1 = 0
                par_y2 = 1

            # we calculate these variables because we will use them to obtain the bilinear interpolation
            Dx, dx = real_partition_x2 - real_partition_x1, par_x - real_partition_x1
            Dy, dy = real_partition_y2 - real_partition_y1, par_y - real_partition_y1

            # we make sure that we do not divide by zero for any case
            if dx == 0:
                Dx = 1
            if dy == 0:
                Dy = 1

            # if par_x2 is bigger than the number of partitions of the pixels of the stamps,
            # make the star be at the first partition of the next pixel
            if par_x2 >= len(self.__stars):
                par_x2 = 0
                pix_x2 = pix_x1 + 1
            else:
                pix_x2 = pix_x1

            # idem for par_y2
            if par_y2 >= len(self.__stars[0]):
                par_y2 = 0
                pix_y2 = pix_y1 + 1
            else:
                pix_y2 = pix_y1

            # stamps necessary for the interpolation
            stamp_11 = self.__stars[par_x1][par_y1]
            stamp_21 = self.__stars[par_x2][par_y1]
            stamp_12 = self.__stars[par_x1][par_y2]
            stamp_22 = self.__stars[par_x2][par_y2]

            # **********************************
            # bilinear interpolation:
            # check the formula at https://es.frwiki.wiki/wiki/Interpolation_bilin%C3%A9aire

            # we generate a matrix in which we will save the star (number of flux on each pixel) after
            # the bilinear interpolation to place it in the image
            star = np.zeros((self.__npx_stamp + 1, self.__npy_stamp + 1))

            if pix_x2 == pix_x1 and pix_y2 == pix_y1:

                for i in range(self.__npx_stamp):
                    for j in range(self.__npy_stamp):
                        star[i][j] = (stamp_21[i][j] - stamp_11[i][j]) * dx / Dx + (
                                      stamp_12[i][j] - stamp_11[i][j]) * dy / Dy + (
                                      stamp_22[i][j] + stamp_11[i][j] - stamp_21[i][j] - stamp_12[i][j]) * dx * dy / (Dx * Dy) + (
                                      stamp_11[i][j])

            elif pix_x2 != pix_x1 and pix_y2 != pix_y1:

                for i in range(self.__npx_stamp + 1):
                    for j in range(self.__npy_stamp + 1):

                        if i - 1 < 0 and j - 1 < 0:
                            star[i][j] = (0.0 - stamp_11[i][j]) * dx / Dx + (
                                          0.0 - stamp_11[i][j]) * dy / Dy + (
                                          0.0 + stamp_11[i][j] - 0.0 -
                                          0.0) * dx * dy / (Dx * Dy) + (
                                          stamp_11[i][j])

                        elif i - 1 < 0 and j >= self.__npy_stamp:
                            star[i][j] = (0.0 - 0.0) * dx / Dx + (
                                          stamp_12[i][j - 1] - 0.0) * dy / Dy + (
                                          0.0 + 0.0 - 0.0 -
                                          stamp_12[i][j - 1]) * dx * dy / (Dx * Dy) + (
                                          0.0)

                        elif i - 1 < 0:
                            star[i][j] = (0.0 - stamp_11[i][j]) * dx / Dx + (
                                          stamp_12[i][j - 1] - stamp_11[i][j]) * dy / Dy + (
                                          0.0 + stamp_11[i][j] - 0.0 -
                                          stamp_12[i][j - 1]) * dx * dy / (Dx * Dy) + (
                                          stamp_11[i][j])

                        elif j - 1 < 0 and i >= self.__npx_stamp:
                            star[i][j] = (stamp_21[i - 1][j] - 0.0) * dx / Dx + (
                                          0.0 - 0.0) * dy / Dy + (
                                          0.0 + 0.0 - stamp_21[i - 1][j] -
                                          0.0) * dx * dy / (Dx * Dy) + (
                                          0.0)

                        elif j - 1 < 0:
                            star[i][j] = (stamp_21[i - 1][j] - stamp_11[i][j]) * dx / Dx + (
                                          0.0 - stamp_11[i][j]) * dy / Dy + (
                                          0.0 + stamp_11[i][j] - stamp_21[i - 1][j] -
                                          0.0) * dx * dy / (Dx * Dy) + (
                                          stamp_11[i][j])

                        elif i >= self.__npx_stamp and j >= self.__npy_stamp:
                            star[i][j] = (0.0 - 0.0) * dx / Dx + (
                                          0.0 - 0.0) * dy / Dy + (
                                          stamp_22[i - 1][j - 1] + 0.0 - 0.0 -
                                          0.0) * dx * dy / (Dx * Dy) + (
                                          0.0)

                        elif i >= self.__npx_stamp:
                            star[i][j] = (stamp_21[i - 1][j] - 0.0) * dx / Dx + (
                                          0.0 - 0.0) * dy / Dy + (
                                          stamp_22[i - 1][j - 1] + 0.0 - stamp_21[i - 1][j] -
                                          0.0) * dx * dy / (Dx * Dy) + (
                                          0.0)

                        elif j >= self.__npy_stamp:
                            star[i][j] = (0.0 - 0.0) * dx / Dx + (
                                          stamp_12[i][j - 1] - 0.0) * dy / Dy + (
                                          stamp_22[i - 1][j - 1] + 0.0 - 0.0 -
                                          stamp_12[i][j - 1]) * dx * dy / (Dx * Dy) + (
                                          0.0)

                        else:
                            star[i][j] = (stamp_21[i - 1][j] - stamp_11[i][j]) * dx / Dx + (
                                          stamp_12[i][j - 1] - stamp_11[i][j]) * dy / Dy + (
                                          stamp_22[i - 1][j - 1] + stamp_11[i][j] - stamp_21[i - 1][j] -
                                          stamp_12[i][j - 1]) * dx * dy / (Dx * Dy) + (
                                          stamp_11[i][j])
            elif pix_x2 != pix_x1:

                for i in range(self.__npx_stamp + 1):
                    for j in range(self.__npy_stamp):

                        if i - 1 < 0:
                            star[i][j] = (0.0 - stamp_11[i][j]) * dx / Dx + (
                                          stamp_12[i][j] - stamp_11[i][j]) * dy / Dy + (
                                          0.0 + stamp_11[i][j] - 0.0 -
                                          stamp_12[i][j]) * dx * dy / (Dx * Dy) + (
                                          stamp_11[i][j])

                        elif i >= self.__npx_stamp:
                            star[i][j] = (stamp_21[i - 1][j] - 0.0) * dx / Dx + (
                                          0.0 - 0.0) * dy / Dy + (
                                          stamp_22[i - 1][j] + 0.0 - stamp_21[i - 1][j] -
                                          0.0) * dx * dy / (Dx * Dy) + (
                                          0.0)

                        else:
                            star[i][j] = (stamp_21[i - 1][j] - stamp_11[i][j]) * dx / Dx + (
                                          stamp_12[i][j] - stamp_11[i][j]) * dy / Dy + (
                                          stamp_22[i - 1][j] + stamp_11[i][j] - stamp_21[i - 1][j] -
                                          stamp_12[i][j]) * dx * dy / (Dx * Dy) + (
                                          stamp_11[i][j])

            elif pix_y2 != pix_y1:

                for i in range(self.__npx_stamp):
                    for j in range(self.__npy_stamp + 1):

                        if j - 1 < 0:
                            star[i][j] = (stamp_21[i][j] - stamp_11[i][j]) * dx / Dx + (
                                          0.0 - stamp_11[i][j]) * dy / Dy + (
                                          0.0 + stamp_11[i][j] - stamp_21[i][j] -
                                          0.0) * dx * dy / (Dx * Dy) + (
                                          stamp_11[i][j])

                        elif j >= self.__npy_stamp:
                            star[i][j] = (0.0 - 0.0) * dx / Dx + (
                                          stamp_12[i][j - 1] - 0.0) * dy / Dy + (
                                          stamp_22[i][j - 1] + 0.0 - 0.0 -
                                          stamp_12[i][j - 1]) * dx * dy / (Dx * Dy) + (
                                          0.0)

                        else:
                            star[i][j] = (stamp_21[i][j] - stamp_11[i][j]) * dx / Dx + (
                                          stamp_12[i][j - 1] - stamp_11[i][j]) * dy / Dy + (
                                          stamp_22[i][j - 1] + stamp_11[i][j] - stamp_21[i][j] -
                                          stamp_12[i][j - 1]) * dx * dy / (Dx * Dy) + (
                                          stamp_11[i][j])

            # **********************************

            for i in range(self.__npx_stamp + 1):
                for j in range(self.__npy_stamp + 1):

                    # we place the matrix (containing the star calculated with the interpolation) in the corresponding
                    # spot of the image by adding (doing a sum) the values of the matrix in each pixel
                    # to the corresponding pixels in the image (initially with value zero)

                    # we make sure that we only add flux if we are inside the image
                    if pix_x1 - int((self.__npx_stamp - 1) / 2) + i in range(self.__npx) and \
                            pix_y1 - int((self.__npy_stamp - 1) / 2) + j in range(self.__npy):
                        self.__array[pix_x1 - int((self.__npx_stamp - 1) / 2) + i][pix_y1 - int((self.__npy_stamp - 1) / 2) + j] += flux * star[i][j]

                    else:
                        pass

    def place_randomStar(self, number_stars, mag1, mag2):
        """This function places random stars in the image. The number of placed stars corresponds to 'number_stars'.
        We choose an interval of magnitudes [mag1, mag2] for the stars."""

        # we generate a catalog that we will fill with random positions and magnitudes (between the chosen interval)
        random_catalog = np.array([[0.] * 3 for i in range(number_stars)])

        for i in range(number_stars):
            # we generate random positions and magnitudes for each star
            rand_x = np.random.uniform(0, self.__npx - 1)
            rand_y = np.random.uniform(0, self.__npy - 1)
            rand_magnitude = np.random.uniform(float(mag1), float(mag2))

            # we add the random values to the catalog
            random_catalog[i] = np.array([rand_x, rand_y, rand_magnitude])

        # we place the stars
        self.placeStar(random_catalog)

    def plotImage(self):
        """This function does a plot of the image."""
        logarithmic_array = np.log10(self.__array + 1.0)  # afegim el +1 pq no calculi log(0) = -inf
        plot = plt.matshow(logarithmic_array)  # norm=colors.Normalize(vmin=0.0, vmax=1.0)
        plt.gray()
        plt.colorbar(plot)
        plt.show()

    def SaveToFile(self):
        """This function saves the image in a file with the name 'GeneratedImage.npy'"""
        np.save('GeneratedImage', self.__array)

