import numpy as np


class Function:
    """
        A class used to calculate different functions

        Attributes
        ----------
        func : int
            used to choose between different possible functions
        pars : list
            parameters of the function (default [0, 0, 0.25, 0.25])

        Methods
        -------
        eval(x,y)
            Evaluates the point (x,y) using the selected function
        """

    def __init__(self, func, pars=[0, 0, 0.25, 0.25]):
        """
        Parameters
        ----------
        func : int
            used to choose between different possible functions:
            func=0 -> normalized gaussian, func=1 -> normalized lorentzian, func=2 -> normalized circle
        pars : list
            parameters of the function: pars=[mu x, mu y, ] (default [0, 0, 0.25, 0.25])
        """

        """This class contains different functions that can be selected using 'func'. For 'func'=0 we have a normalized gaussian,
        for 'func'=1 we have a normalized lorentzian, for 'func'=2 we have a normalized circle.
        'pars' is a list used to change the function parameters. The different components are:
        pars (= parameters) components --> 1st & 2nd components:     {mu_x, mu_y};
                                           3rd & 4th components:     {sigma_x, sigma_y} for func = 0;
                                                                     {gamma_x, gamma_y} for func = 1;
                                                                     {radius} for func = 2"""

        self.func = func
        self.pars = pars

        self.__val = 0.0

    def eval(self, x, y):
        """Evaluate the point (x,y) with the chosen function (func = 0, 1, 2)"""

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
