import numpy as np
import matplotlib.pyplot as plt
from StampCreator import StampCreator
from Stamp import Stamp
import os
import sys


class StampLibrary:
    """
    A class used to generate a library of stamps placing a star in different positions of the central pixel

    This class uses the modules StampCreator, Stamp, numpy, matplotlib.pyplot, os and sys

    Attributes
    ----------
    func: int
            used to choose between different possible functions (from Functions module):
            func=0 -> normalized gaussian, func=1 -> normalized lorentzian, func=2 -> normalized circle
            (see Function module documentation for more information)
    npar_x : int
            number of partitions of the central pixel of the stamp on the x side (we move the star through each partition)
    npar_y : int
            number of partitions of the central pixel of the stamp on the y side (we move the star through each partition)
    npix_x: int
            number of pixels of the stamp on the x-axis
    npix_y: int
            number of pixels of the stamp on the y-axis

    Methods
    -------
    generateLibrary()
        Generates the library of stamps

    SaveToFile()
        Saves the library of stamps to a file in .npy format

    DrawLibrary()
        Does a multiplot of all the stamps of the library of stamps

    getError()
        Obtains the maximum difference of flux inside a stamp with respect to flux = 1.0
        (supposing that the flux inside all the stamps is normalized to 1.0)
    """

    def __init__(self, func, npar_x, npar_y, npix_x, npix_y):
        """Initialize the class StampLibrary

        :param int func: used to choose between different possible functions (from Functions module):
                func=0 -> normalized gaussian, func=1 -> normalized lorentzian, func=2 -> normalized circle
                (see Function module documentation for more information)
        :param int npar_x : number of partitions of the central pixel of the stamp on the x side (we move the star through each partition)
        :param int npar_y : number of partitions of the central pixel of the stamp on the y side (we move the star through each partition)
        :param int npix_x: number of pixels of the stamp on the x-axis
        :param int npix_y: number of pixels of the stamp on the y-axis
        """

        self.func = func

        # number of partitions of the central pixel
        self.par_x = npar_x
        self.par_y = npar_y

        # number of pixels of the stamp
        self.pix_x = npix_x
        self.pix_y = npix_y

        # initialize the creator of stamps with the chosen function and number of pixels
        self.generated_stamp = StampCreator(self.func, self.pix_x, self.pix_y)

        # calculate longitude of the partitions of the central pixel
        self.x_shift = self.generated_stamp.lx / self.par_x
        self.y_shift = self.generated_stamp.ly / self.par_y

        # create an empty library that will be filled with stamps
        self.Library = [[0.]*self.par_x for i in range(self.par_y)]

    def generateLibrary(self):
        """Generate the library of stamps (placing a unique star in all the partitions of the central pixel)

        :return: Saves the generated library in a class variable
        """

        print('\n This process can take several minutes... \n')

        for i in range(self.par_x):
            for j in range(self.par_y):

                # we generate a new stamp each time the star is placed in a new partition
                LibraryElement = self.generated_stamp.generateStamp(-self.generated_stamp.lx / 2.0 + self.x_shift / 2.0 + self.x_shift * i, -self.generated_stamp.ly / 2.0 + self.y_shift / 2.0 + self.y_shift * j)

                # we save each stamp in a library (class variable)
                self.Library[i][j] = LibraryElement

    def SaveToFile(self):
        """Save the library of stamps to a file in .npy format in a folder called 'stamp_libraries'

        :return: Saves the .npy file in a folder called 'stamp_libraries' with the name p.ex. 'Library_5x5.ny'
        in the case we have stamps of 5x5 pixels
        """

        # path were the library is saved
        output_path = sys.path[0] + '/stamp_libraries'
        # name of the library file
        file_name = 'Library_'+str(self.pix_x)+"x"+str(self.pix_y)
        # the library is saved
        np.save(os.path.join(output_path, file_name), self.Library)

    def DrawLibrary(self):
        """Do a multiplot of all the stamps of the library of stamps

        :return: A plot of all the stamps of the library together using matplotlib.pyplot.subplots
        and matplotlib.pyplot.imshow
        """

        # call a function to do multiplots
        fig, axs = plt.subplots(self.par_x, self.par_y)
        # set the title of the plot
        fig.suptitle('Library')

        # move through the stamps of the library (x-axis)
        for i in range(self.par_x):
            # move through the stamps of the library (y-axis)
            for j in range(self.par_y):
                # for each step of the loop plot the corresponding stamp of the library
                axs[i, j].imshow(self.Library[i][j])
                # place label axis only on the outer stamps of the multiplot
                axs[i, j].label_outer()

        # plot the multiplot
        plt.show()

    def getError(self):
        """Obtain the maximum difference of flux inside a stamp with respect to flux = 1.0
        (supposing that the flux inside all the stamps is normalized to 1.0)

        :return: the maximum error of flux (the one that is farther from 1.0)
        """

        # we create an array were we will save the value of flux of all the stamps of the library
        finalSum = np.zeros([self.par_x, self.par_y])

        # we move through all the elements of the library (all the stamps)
        for i in range(self.par_x):
            for j in range(self.par_y):

                # we do the sum of all the flux inside each stamp and save the value of each of them as a new
                # element on the array we created (we use the module Stamp.py)
                finalSum[i][j] = Stamp(self.Library[i][j]).getSum()

        # we obtain the maximum and minimum value of flux of all the stamps of the library
        # (the one stamp with the lower value and the one with the higher)
        maxValue = np.max(finalSum)
        minValue = np.min(finalSum)

        # we calculate the difference of the min and max value with respect to 1.0
        # (1.0 is the value of flux that the stamp should have because it is normalized)
        maxError = abs(1.0-maxValue)
        minError = abs(1.0-minValue)

        # if the error of the max value is greater than the error of the min value, the final error is the
        # error of the max value
        if maxError >= minError:
            Error = maxError
        # if the error of the min value is greater than the error of the max value, the final error is the
        # error of the min value
        elif minError > maxError:
            Error = minError

        # we print the final error
        print('The maximum error is: ', Error)
