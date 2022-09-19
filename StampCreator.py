import numpy as np
from Function import Function


class StampCreator:
    """
    A class used to create stamps (matrix of pixels)

    This class uses the modules Function and numpy

    Attributes
    ----------
    func: int
            used to choose between different possible functions (from Functions module):
            func=0 -> normalized gaussian, func=1 -> normalized lorentzian, func=2 -> normalized circle
            (see Function module documentation for more information)
    px : int
            number of pixels of the stamp in the x side
    py : int
            number of pixels of the stamp in the y side

    Methods
    -------
    setNsubpix(nx_input, ny_input)
        Sets the number of partitions inside each pixel, this net used to calculate the value of the chosen function on each point

    generateStamp(x_shift, y_shift)
        Generates a stamp using the chosen function (with 'func') and the desired number of pixels (with 'px' & 'py')
    """

    # partitions inside each pixel
    # (used to set a smaller net inside the pixel were we calculate the value of the chosen function on each point)
    __nx = 100  # partitions on x
    __ny = 100  # partitions on y

    # longitude of the pixels (set to unity)
    lx = 1.0  # long. pixel on x
    ly = 1.0  # long. pixel on y
    # note: those are not protected variables bc are needed in another module

    # longitude of the partitions of the net inside the pixel
    __dx = lx / __nx
    __dy = ly / __ny

    def __init__(self, func, px, py):
        """Initialize the class StampCreator

        :param int func: used to choose between different possible functions (from Functions module):
                func=0 -> normalized gaussian, func=1 -> normalized lorentzian, func=2 -> normalized circle
                (see Function module documentation for more information)
        :param int px : number of pixels of the stamp in the x side
        :param int py : number of pixels of the stamp in the y side
        """

        self.func = func

        self.px = px
        self.py = py

    def setNsubpix(self, nx_input, ny_input):
        """Set the number of partitions inside each pixel (this net is used to calculate the value of the chosen function on each point)

        :param int nx_input: partitions inside the pixel in the x side (default 100)
        :param int ny_input: partitions inside the pixel in the y side (default 100)
        """

        self.__nx = nx_input
        self.__ny = ny_input

    def generateStamp(self, x_shift=0.0, y_shift=0.0):
        """Generate a stamp using the chosen function (with 'func') and the desired number of pixels (with 'px' & 'py')

        :param float x_shift: shift of the origin of the chosen function on the x-axis (default 0.0)
        :param float y_shift: shift of the origin of the chosen function on the y-axis (default 0.0)

        :return: generated stamp
        :rtype: array
        """

        # initialize the array that will contain the values of the stamp on each pixel
        myStamp = np.zeros([self.px, self.py])

        # move through the x-axis of the pixels of the stamp
        for I in range(self.px):
            # move through the y-axis of the pixels of the stamp
            for J in range(self.py):
                # initialize the value inside each pixel of the stamp (and set to zero after each loop)
                val = 0.0
                # move through the x-axis of the net inside the pixel
                for i in range(self.__nx):
                    # move through the y-axis of the net inside the pixel
                    for j in range(self.__ny):
                        # calculate the x & y positions each time that we are moving through the net

                        # x = (coords. origin) + (displacement of 1 pixel on each iteration) +
                        # + (we place the x in the middle of the pixel) + (we move through the small net
                        # on each iteration) + (we add a shift from the origin)

                        x = - (self.px * self.lx) / 2.0 + I * self.lx + self.__dx / 2.0 + i * self.__dx - x_shift
                        y = - (self.py * self.ly) / 2.0 + J * self.ly + self.__dy / 2.0 + j * self.__dy - y_shift

                        # evaluate the function (of Function module) on each point of the net and sum all results
                        # inside 'val' (we only save the final value of the function inside the pixel, not the smaller net)
                        val += Function(self.func).eval(x, y)

                # save the value of each pixel in myStamp (array of pixels)
                # the final value is normalized (multiplying by self.__dx & self.__dy)
                myStamp[I, J] = val * self.__dx * self.__dy

        return myStamp

