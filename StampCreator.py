import numpy as np
from Function import Function


class StampCreator:
    """
    A class used to create stamps (matrix of pixels)

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

    generateStamp(x_shift=0.0, y_shift=0.0)
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

        Parameters
        ----------
        func: int
                used to choose between different possible functions (from Functions module):
                func=0 -> normalized gaussian, func=1 -> normalized lorentzian, func=2 -> normalized circle
                (see Function module documentation for more information)
        px : int
                number of pixels of the stamp in the x side
        py : int
                number of pixels of the stamp in the y side
        """

        self.func = func

        self.px = px
        self.py = py

    def setNsubpix(self, nx_input, ny_input):
        """Set the number of partitions inside each pixel (this net is used to calculate the value of the chosen function on each point)

        Parameters
        ----------
        nx_input: int
                partitions inside the pixel in the x side (default 100)
        ny_input: int
                partitions inside the pixel in the y side (default 100)
        """

        self.__nx = nx_input
        self.__ny = ny_input

    def generateStamp(self, x_shift=0.0, y_shift=0.0):
        """Generate a stamp using the chosen function (with 'func') and the desired number of pixels (with 'px' & 'py')

        Parameters
        ----------
        x_shift: float
                shift of the origin of the chosen function on the x-axis
        y_shift: float
                shift of the origin of the chosen function on the y-axis

        Returns
        -------
        myStamp: array
            generated stamp
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
                        x = - (self.px * self.lx) / 2.0 + I * self.lx + self.__dx / 2.0 + i * self.__dx - x_shift  # posem (-) x_shift pq així, al posar un desplaçament positiu, la funció es deslpaça cap a un numero mes gran
                        y = - (self.py * self.ly) / 2.0 + J * self.ly + self.__dy / 2.0 + j * self.__dy - y_shift

                        # evaluate the function (of Function module) on each point of the net and sum all results
                        # inside 'val' (we only save the final value of the function inside the pixel, not the smaller net)
                        val += Function(self.func).eval(x, y)

                # save the value of each pixel in myStamp (array of pixels)
                # the final value is normalized (multiplying by self.__dx & self.__dy)
                myStamp[I, J] = val * self.__dx * self.__dy

        return myStamp

