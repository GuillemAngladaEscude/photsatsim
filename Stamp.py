# import numpy as np
import matplotlib.pyplot as plt
# from StampCreator import StampCreator as sc


class Stamp:

    __scalefactor = 1.0

    def __init__(self, myStamp):

        # number of pixels
        self.px = len(myStamp)
        self.py = len(myStamp[0])

        # introduced stamp
        self.myStamp = myStamp

    def getSum(self):
        sum = 0.0

        for i in range(self.px):
            for j in range(self.py):

                sum += self.myStamp[i, j] * self.__scalefactor

        return sum

    def setScale(self, fac):
        self.__scalefactor = fac

    def drawStamp(self):
        plt.matshow(self.myStamp)
        plt.show()


'''my_stamp = sc(5,5).generateStamp(0)
# my_stamp = np.array([[0,0,0,0,0], [1,2,3,4,5], [5,4,3,2,1], [3,3,3,3,3], [0,0,1,1,0], [0,0,0,0,0]])

print(Stamp(my_stamp).getSum())
do_something = Stamp(my_stamp).drawStamp()'''

