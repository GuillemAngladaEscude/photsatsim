import numpy as np
import matplotlib.pyplot as plt
import StampCreator as sc
import Stamp as s
import os
import sys


class StampLibrary:

    def __init__(self, npar_x, npar_y, npix_x, npix_y):
        self.par_x = npar_x   # number of partitions on x
        self.par_y = npar_y

        self.pix_x = npix_x   # number of pixels in the stamp
        self.pix_y = npix_y

        self.stamp = sc.StampCreator(self.pix_x, self.pix_y)

        self.x_shift = self.stamp.lx / self.par_x
        self.y_shift = self.stamp.ly / self.par_y

        self.Library = [[0.]*self.par_x for i in range(self.par_y)]

    def generateLibrary(self):

        for i in range(self.par_x):
            for j in range(self.par_y):

                LibraryElement = self.stamp.generateStamp(0, -self.stamp.lx / 2.0 + self.x_shift / 2.0 + self.x_shift * i, -self.stamp.ly / 2.0 + self.y_shift / 2.0 + self.y_shift * j)
                self.Library[i][j] = LibraryElement

    def SaveToFile(self):
        output_path = sys.path[0] + '/stamp_libraries'
        file_name = 'Library_'+str(self.pix_x)+"x"+str(self.pix_y)
        np.save(os.path.join(output_path, file_name), self.Library)

    def DrawLibrary(self):
        # Library_ = np.load('Library_'+str(self.pix_x)+"x"+str(self.pix_y)+'.npy')  # todo: problema -> primer has de guardar el fitxer perq sino estaras imprimint una cosa que no te pq coincidir amb el que li demanes

        fig, axs = plt.subplots(self.par_x, self.par_y)
        fig.suptitle('Library')

        # images = [[0.] * self.par_x for i in range(self.par_y)]

        for i in range(self.par_x):
            for j in range(self.par_y):
                # images = axs[i, j].imshow(self.Library[i][j])
                axs[i, j].imshow(self.Library[i][j])
                axs[i, j].label_outer()

        plt.show()

    def getError(self):  # todo: pensar si cridem la llibreria aixi. pot ser que volguem cridar una llibreria pel seu nom.
        # Library_ = np.load('Library_'+str(self.pix_x)+"x"+str(self.pix_y)+'.npy')

        finalSum = np.zeros([self.par_x, self.par_y])

        for i in range(self.par_x):
            for j in range(self.par_y):

                finalSum[i][j] = s.Stamp(self.Library[i][j]).getSum()

        maxValue = np.max(finalSum)
        minValue = np.min(finalSum)

        maxError = abs(1.0-maxValue)
        minError = abs(1.0-minValue)

        # print(maxValue, minValue, maxError, minError)

        if maxError >= minError:
            Error = maxError
        elif minError > maxError:
            Error = minError

        print('The maximum error is: ', Error)













