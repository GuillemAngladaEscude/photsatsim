import numpy as np
import matplotlib.pyplot as plt
import sys
import os


class Image:
    """
        A class used to generate an image using a library of stamps to place the stars

        If 'stampSize' is introduced, a library of stamps with that number of pixels on each side is selected
            (the default is stamps of 5 pixels on each side).
        If 'libname' is introduced, the image will be generated using the library of stamps that has that name
            (this needs only to be used if we want to work with a library that is not generated with 'Driver_library.py')

        This class uses the modules numpy, matplotlib.pyplot, os and sys

        Attributes
        ----------
        npx: int
                number of pixels of the image on the x-axis
        npy: int
                number of pixels of the image on the y-axis
        stampSize: int (default 5)
                number of pixels of the stamps (supposed to have the same number of pixels on both sides)
        libname: srt (default 'none')
                name of the library of stamps that we want to work with

        Methods
        -------
        placeStar(catalog)
            Places each star of the catalog on the image

        SaveToFile()
            Saves the generated image to a file in .npy format

        plotImage()
            Does a plot of the generated image in logarithmic scale

        place_randomStar(number_stars, mag1, mag2)
            Places stars in random positions and with random apparent magnitudes (inside an interval) on the image

        get_array()
            Returns the generated image saved in memory (which is a private class variable)

        get_arrayCopy()
            Generates a copy of the generated image

        test(number_partitions, pixel_x_input, pixel_y_input)
            Tests that the flux inside the image remains constant by placing a unique star in different points of a pixel
        """

    def __init__(self, npx, npy, stampSize = 5, libname='none'):
        """Initialize the class Image

        :param int npx: number of pixels of the image on the x-axis
        :param int npy: number of pixels of the image on the y-axis
        :param int stampSize: odd int that sets the number of pixels of the stamps
        (supposed to have the same number of pixels on both sides) (default 5)
        :param srt libname: name of the library of stamps that we want to work with (default 'none')
        """

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

    def __setPixelsStamp(self, npx_stamp_input, npy_stamp_input):
        """Sets the number of pixels that the stamps have on each side

        :raises ValueError: if the stamp does not have an odd number of pixels on each side
        """

        self.__npx_stamp = npx_stamp_input
        self.__npy_stamp = npy_stamp_input
        # The stamp must have an odd number of pixels at each side. If this does not happen, an error occurs
        if self.__npx_stamp % 2 == 0 or self.__npy_stamp % 2 == 0:
            raise ValueError('The stamp must have an odd number of pixels at each side.')

    def __loadLibrary(self, libraryname):
        """Loads a file corresponding to a library of stamps with the name 'libraryname' stored in the folder 'stamp_libraries'.
        This library is used to generate the image
        """

        output_path = sys.path[0] + '/stamp_libraries/'
        self.__stars = np.load(os.path.join(output_path, libraryname))

    def __initLibraryWithStampSize(self, stampSize):
        """Initializes a library with the stamps of the size of
        'stampSize' (supposed squares), which shall be an odd integer (1,3,5,7,9...)
        """

        libname = 'Library_'+str(stampSize)+"x"+str(stampSize)+'.npy'
        self.__setPixelsStamp(stampSize, stampSize)
        self.__loadLibrary(libname)

    def __initLibraryWithLibraryName(self, libname_input):
        """Initializes a library with the name introduced in 'libname_input'.
        Then obtains the number of pixels of the stamps of the library (supposed squares).
        Finally, sets the number of pixels of the stamps
        """

        self.__loadLibrary(libname_input)
        stampSize = self.__stars[0][0]
        self.__setPixelsStamp(stampSize, stampSize)

    def __mag_to_flux(self, magnitude):
        """Receives an apparent magnitude and returns the corresponding
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
        """Return the generated image (which is a private class variable). Modifying this array means modify the image
        saved in memory

        :return: The private class variable 'self.__array' which is the generated image
        """

        return self.__array

    def get_arrayCopy(self):
        """Create a copy of the generated image. Modifying this array has no effect on the original image saved in memory

        :return: A copy of the image saved in memory 'self.__array'
        """

        copy = self.__array.copy()
        return copy

    def test(self, number_partitions, pixel_x_input=4.378, pixel_y_input=5):
        """Test that the flux inside the image remains constant while we place a unique star
        (with apparent magnitude = 5.0) in different points of a pixel

        The x-axis is fixed and a star is moved through the y-axis from one end of the pixel to the other.

        :param number_partitions: number of partitions that we want to do of the pixel for the test

        :param pixel_x_input: float that sets the x-position of the star on the image (fixes the x-axis) (default 4.378)

        :param pixel_y_input: int that sets the pixel on the y-axis were we want to move the star (default 5)

        :return: The main difference between the maximum and minimum value of flux and a plot of all the values of flux
        on each partition of the pixel were a star is placed
        """

        # position un x and y were we want to place the star (x-axis fixed and y-axis the star moves)
        pixel_x, pixel_y = pixel_x_input, pixel_y_input

        # length of the partitions of the pixel on the y-axis
        dy = 1.0 / number_partitions

        # we create the arrays that will contain the value of the flux on each partition (results)
        # and the one that will contain the positions it corresponds on the pixel (positions)
        results = np.zeros(number_partitions)
        positions = np.zeros(number_partitions)

        # loop that moves through the y-axis in each partition of the pixel
        for j in range(number_partitions):

            # we create the variable that will contain the sum of the flux inside the whole image
            # (it is set to zero each time the star is placed in the next partition to do the sum again)
            flux_sum = 0.0

            # position of the star in the image (its pixel position + partition) in the y-axis
            position_y = pixel_y + j * dy

            # we create a catalog that we use to place the star each time
            test_catalog = [[pixel_x, position_y, 5.0]]

            # we save each new position on the y-axis of the pixel on a new element of the list that is created earlier
            positions[j] = position_y

            # we generate a new image each time we move the star
            test_image = Image(self.__npx, self.__npy)
            # we place the star in the corresponding spot (with 'test_catalog')
            test_image.placeStar(test_catalog)

            # we move through the pixels of the image
            for I in range(self.__npx):
                for J in range(self.__npy):
                    # we sum all the flux on each pixel of the image and obtain the total flux on the image
                    flux_sum += test_image.__array[I][J]

            # we save the total flux of the image as a new element of the list created to save the results
            # (each image saves its flux as a new element of the list)
            results[j] = flux_sum

        # once we have all the results, we obtain the maximum and minimum values
        maximum = max(results)
        minimum = min(results)

        # we calculate the difference between the maximum and minimum value
        main_difference = maximum - minimum
        # print the main difference
        print('The main difference between different positions of the star is:', main_difference)

        # do a plot of all the values of flux with respect of the positions of the stars
        plt.plot(positions, results)
        plt.xlabel('pixel position')
        plt.ylabel('flux')
        plt.show()

    def placeStar(self, catalog):
        """Place each star of 'catalog' on the image

        The stars are placed doing a bilinear interpolation of the stamps (see: https://es.frwiki.wiki/wiki/Interpolation_bilin%C3%A9aire)

        :param list catalog: contains the position of the stars and their apparent magnitude with the following shape
                        [[position pixel on x of star 1, position pixel on y of star 1, apparent magnitude of star 1], [idem for star 2], [idem for star 3], ...]
                        (list of lists)

        :return: saves the generated image in memory (is a class variable)
        """

        for k in range(len(catalog)):

            # we assign each value of the catalog to its corresponding variable
            x, y, flux = catalog[k][0], catalog[k][1], self.__mag_to_flux(catalog[k][2])

            x, y = self.xy_2_xy_TL(x,y)  # afegit per Lluc

            # we calculate the pixel and the partition of the pixel in which is the star
            pix_x1, par_x = int(x // 1.0), x % 1.0
            pix_y1, par_y = int(y // 1.0), y % 1.0

            # we calculate the number of partitions that have the pixels of the stamps (the precision of the library)
            partitions_x_stamp = len(self.__stars)
            partitions_y_stamp = len(self.__stars[0])

            # with the number of partitions we calculate the length of each partition
            dx_partitions = 1.0 / float(partitions_x_stamp)
            dy_partitions = 1.0 / float(partitions_y_stamp)

            # with this loop we identify between which stamps is our star (for the bilinear interpolation) in the x-axis
            for i in range(partitions_x_stamp + 1):
                if par_x <= (i * dx_partitions):
                    break

            # we assign the result of the loop to these variables in order to use it later to call the correct stamps
            par_x1 = i - 1
            par_x2 = i

            # we calculate these variables because are necessary to compute the bilinear interpolation
            real_partition_x1 = par_x1 * dx_partitions
            real_partition_x2 = par_x2 * dx_partitions

            # now we make sure that if the star is in the 0.0 spot, the stamp is correctly assigned
            if par_x1 == -1:
                par_x1 = 0
                par_x2 = 1

            # with this loop we identify between which stamps is our star (for the bilinear interpolation) in the y-axis
            for j in range(partitions_y_stamp + 1):
                if par_y <= (j * dy_partitions):
                    break

            # we assign the result of the loop to these variables in order to use it later to call the correct stamps
            par_y1 = j - 1
            par_y2 = j

            # we calculate these variables because are necessary to compute the bilinear interpolation
            real_partition_y1 = par_y1 * dy_partitions
            real_partition_y2 = par_y2 * dy_partitions

            # now we make sure that if the star is in the 0.0 spot, the stamp is correctly assigned
            if par_y1 == -1:
                par_y1 = 0
                par_y2 = 1

            # we calculate these variables because we will use them to obtain the bilinear interpolation
            Dx, dx = real_partition_x2 - real_partition_x1, par_x - real_partition_x1
            Dy, dy = real_partition_y2 - real_partition_y1, par_y - real_partition_y1

            # we make sure that we do not divide by zero for any case (see the formula of the link)
            if dx == 0:
                Dx = 1
            if dy == 0:
                Dy = 1

            # if par_x2 is bigger than the number of partitions of the pixels of the stamps (it is outside the pixel),
            # make the star be at the first partition of the next pixel
            if par_x2 >= len(self.__stars):
                par_x2 = 0
                pix_x2 = pix_x1 + 1
            # if it is not bigger, both stamps are on the same pixel
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

            # if the star is between stamps that are on the same pixel (for x and y) do the following
            if pix_x2 == pix_x1 and pix_y2 == pix_y1:

                for i in range(self.__npx_stamp):
                    for j in range(self.__npy_stamp):
                        # apply the bilinear interpolation formula (see the link)
                        star[i][j] = (stamp_21[i][j] - stamp_11[i][j]) * dx / Dx + (
                                      stamp_12[i][j] - stamp_11[i][j]) * dy / Dy + (
                                      stamp_22[i][j] + stamp_11[i][j] - stamp_21[i][j] - stamp_12[i][j]) * dx * dy / (Dx * Dy) + (
                                      stamp_11[i][j])

            # if the star is between stamps that are on different pixels -one pixel and the next one- (for x and y) do the following
            elif pix_x2 != pix_x1 and pix_y2 != pix_y1:

                for i in range(self.__npx_stamp + 1):
                    for j in range(self.__npy_stamp + 1):

                        # apply the bilinear interpolation formula taking zero value for all the elements that are
                        # outside the stamps (we should call a non existing matrix element, as a solution we say
                        # that the value of the element is zero)

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

            # if the star is between stamps that are on different pixels -one pixel and the next one- (only for x) do the following
            elif pix_x2 != pix_x1:

                for i in range(self.__npx_stamp + 1):
                    for j in range(self.__npy_stamp):

                        # apply the bilinear interpolation formula taking zero value for all the elements that are
                        # outside the stamps (we should call a non existing matrix element, as a solution we say
                        # that the value of the element is zero)

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

            # if the star is between stamps that are on different pixels -one pixel and the next one- (only for y) do the following
            elif pix_y2 != pix_y1:

                for i in range(self.__npx_stamp):
                    for j in range(self.__npy_stamp + 1):

                        # apply the bilinear interpolation formula taking zero value for all the elements that are
                        # outside the stamps (we should call a non existing matrix element, as a solution we say
                        # that the value of the element is zero)

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
        """Place stars in random positions and with random apparent magnitudes (inside an interval) on the image

        :param number_stars: number of placed stars
        :param mag1: first element of the interval of magnitudes
        :param mag2: last element of the interval of magnitudes

        :return: saves the generated image in memory (is a class variable)
        """

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
        """Do a plot of the image in logarithmic scale (base 10)

        :return: a plot of the generated image using matplotlib.pyplot.matshow
        """

        # we do a change of scale from linear to logarithmic (base 10)
        # we add a + 1.0 to the generated image to make sure that we do not calculate log(0) = -inf
        logarithmic_array = np.log10(self.__array + 1.0)

        # we do a plot of the logarithmic image
        plot = plt.matshow(logarithmic_array)
        # colors of the plot are gray
        plt.gray()

        # the plot has a color bar
        plt.colorbar(plot)

        plt.show()

    def SaveToFile(self, i):
        """Save the image in a file inside the folder 'output'

        :return: the generated image is saved in a file called 'GeneratedImage.npy' inside the folder 'output'
        """
        logarithmic_array = np.log10(self.__array + 1.0)

        # we do a plot of the logarithmic image

        plt.figure(figsize=(10, 10))

        plt.matshow(logarithmic_array, fignum=1)
        plt.axis('off')

        # colors of the plot are gray
        plt.gray()

        #plot(figsize=(8, 8), dpi=80)
        #plt.show()

        # we make sure that the image is saved to a folder called 'output' in the same directory
        output_path = sys.path[0] + '/output/'
        file_name = 'FullScan_Image' + str(i)

        #np.save(os.path.join(output_path, file_name), self.__array)
        plt.savefig(output_path + file_name, bbox_inches="tight")

        plt.close()


    def xy_2_xy_TL(self, x, y):
        return [self.__npx/2 - y, self.__npy/2 + x]

