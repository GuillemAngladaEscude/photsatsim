import numpy as np


class Function:
    """
    A class used to calculate a point (x,y) with different possible functions

    This class uses the module numpy

    Attributes
    ----------
    func : int
            used to choose between different possible functions:
            func=0 -> normalized gaussian, func=1 -> normalized lorentzian, func=2 -> normalized circle
    pars : list
            parameters of the function (default [0, 0, 0.25, 0.25])
            normalized gaussian: pars=[mu_x, mu_y, sigma_x, sigma_y]
            normalized lorentzian: pars=[mu_x, mu_y, gamma_x, gamma_y]
            normalized circle: pars=[mu_x, mu_y, radius]

    Methods
    -------
    eval(x,y)
        Evaluates the point (x,y) using the selected function (with 'func')
    """

    def __init__(self, func, pars=[0, 0, 0.25, 0.25]):
        """Initialize the class Function

        :param int func : used to choose between different possible functions:
        func=0 -> normalized gaussian, func=1 -> normalized lorentzian, func=2 -> normalized circle

        :param list pars : parameters of the function (default [0, 0, 0.25, 0.25]). For a normalized gaussian:
        pars=[mu_x, mu_y, sigma_x, sigma_y]. For a normalized lorentzian: pars=[mu_x, mu_y, gamma_x, gamma_y].
        For a normalized circle: pars=[mu_x, mu_y, radius]
        """

        self.func = func
        self.pars = pars

        # final value of the function after evaluating on the chosen point (x,y)
        # initialized with value = 0.0
        self.__val = 0.0

    def eval(self, x, y):
        """Evaluate the point (x,y) with the chosen function

        :param float x: x position in which we want to calculate the function
        :param float y: y position in which we want to calculate the function

        :return: the value of the function in the point (x, y)
        :rtype: float
        """

        # definition of the variable 'val', that gives the value of the function in each point (x, y)
        self.__val = 0.0

        # func == 0 --> normalized gaussian:
        if self.func == 0:
            self.__val = (1.0/(self.pars[2]*self.pars[3]*2*np.pi))*np.exp(-((((x-self.pars[0])**2)/(2*self.pars[2]**2))+(((y-self.pars[1])**2)/(2*self.pars[3]**2))))

        # func == 1 --> normalized lorentzian:
        if self.func == 1:
            self.__val = 1.0/(np.pi**2*self.pars[2]*self.pars[3]*(1.0+((x-self.pars[0])/self.pars[2])**2)*(1.0+((y-self.pars[1])/self.pars[3])**2))

        # func == 2 --> normalized circle:
        if self.func == 2:
            if (x-self.pars[0])**2+(y-self.pars[1])**2 <= (self.pars[2])**2:
                self.__val = 1.0 / (np.pi * self.pars[2]**2)
            else:
                self.__val = 0

        # we obtain the value of the function in a concrete point (x, y)
        return self.__val
