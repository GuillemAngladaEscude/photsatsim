import matplotlib.pyplot as plt


class Stamp:
    """
    A class used to work with a stamp (matrix of pixels)

    Attributes
    ----------
    myStamp : array
            matrix of pixels that we want to work with

    Methods
    -------
    setScale(fac)
        Sets a factor for which the stamp will be multiplied

    getSum()
        Gives the total sum of the values inside each pixel of the stamp

    drawStamp()
        Does a plot of the introduced stamp
    """

    # initialize a scale factor (default 1.0) for which the stamp is multiplied
    # method setScale(fac) can change its value
    __scalefactor = 1.0

    def __init__(self, myStamp):
        """Initialize the class Stamp

        Parameters
        ----------
        myStamp : array
                matrix of pixels that we want to work with
        """

        # number of pixels on each side of the stamp
        self.__px = len(myStamp)  # pixels on x
        self.__py = len(myStamp[0])  # pixels on y

        # introduced stamp
        self.myStamp = myStamp

    def setScale(self, fac):
        """Set the factor for which the stamp will be multiplied

        Parameters
        ----------
        fac : float
                factor for which the stamp is multiplied
        """

        # change the value of the scale factor for the introduced value
        self.__scalefactor = fac

    def getSum(self):
        """Give the total sum of the values inside each pixel of the stamp

        Returns
        -------
        sum: float
                total sum of the values inside each pixel of the stamp
        """

        # initialize the value of the sum (at first = 0.0)
        sum = 0.0

        # move through the x-axis of the stamp (matrix of pixels)
        for i in range(self.__px):
            # move through the y-axis of the stamp
            for j in range(self.__py):

                # sum the value of each pixel of the stamp
                sum += self.myStamp[i, j] * self.__scalefactor

        return sum

    def drawStamp(self):
        """Do a plot of the introduced stamp

        Returns
        -------
        a plot of the introduced stamp using matplotlib.pyplot.matshow
        """

        plt.matshow(self.myStamp)
        plt.show()
