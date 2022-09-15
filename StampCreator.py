import numpy as np
import Function as f


class StampCreator:

    nx = 100  # partitions inside each pixel  #todo: variable?
    ny = 100
    lx = 1.0  # long. pixel
    ly = 1.0
    dx = lx / nx
    dy = ly / ny

    def __init__(self, px, py):
        self.px = px  # num. pixels
        self.py = py

    def setNsubpix(self, nx_input, ny_input):
        self.nx = nx_input
        self.ny = ny_input

    def generateStamp(self, func, x_shift=0.0, y_shift=0.0):
        myStamp = np.zeros([self.px, self.py])

        for I in range(self.px):
            for J in range(self.py):
                val = 0.0
                for i in range(self.nx):
                    for j in range(self.ny):
                        x = - (self.px * self.lx) / 2.0 + I * self.lx + self.dx / 2.0 + i * self.dx - x_shift  # posem (-) x_shift pq així, al posar un desplaçament positiu, la funció es deslpaça cap a un numero mes gran
                        y = - (self.py * self.ly) / 2.0 + J * self.ly + self.dy / 2.0 + j * self.dy - y_shift

                        val += f.Function(func).eval(x, y)

                myStamp[I, J] = val * self.dx * self.dy

        return myStamp

