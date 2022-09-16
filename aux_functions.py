import numpy as np
from numpy import pi, ndarray
from scipy.optimize import root

def R3(theta):  # todo: revisar les tres matrius. Al curtis surten les transposades com si no ho fossin
    R = np.zeros((3, 3))
    R[0, 0] = np.cos(theta)
    R[0, 1] = -np.sin(theta)
    R[1, 0] = np.sin(theta)
    R[1, 1] = np.cos(theta)
    R[2, 2] = 1
    return R


def R2(theta):
    R = np.zeros((3, 3))
    R[0, 0] = np.cos(theta)
    R[0, 2] = np.sin(theta)
    R[1, 1] = 1
    R[2, 0] = -np.sin(theta)
    R[2, 2] = np.cos(theta)
    return R


def R1(theta):
    R = np.zeros((3, 3))
    R[0, 0] = 1
    R[1, 1] = np.cos(theta)
    R[1, 2] = -np.sin(theta)
    R[2, 1] = np.sin(theta)
    R[2, 2] = np.cos(theta)
    return R


def Rs(s, alpha):
    s = s / np.linalg.norm(s)
    sx = s[0]
    sy = s[1]
    sz = s[2]
    Q = np.zeros((3, 3))
    Q[0, 0] = np.cos(alpha) + (1 - np.cos(alpha)) * sx ** 2
    Q[0, 1] = (1 - np.cos(alpha)) * sx * sy - sz * np.sin(alpha)
    Q[0, 2] = (1 - np.cos(alpha)) * sx * sz + sy * np.sin(alpha)
    Q[1, 0] = (1 - np.cos(alpha)) * sx * sy + sz * np.sin(alpha)
    Q[1, 1] = np.cos(alpha) + (1 - np.cos(alpha)) * sy ** 2
    Q[1, 2] = (1 - np.cos(alpha)) * sy * sz - sx * np.sin(alpha)
    Q[2, 0] = (1 - np.cos(alpha)) * sx * sz - sy * np.sin(alpha)
    Q[2, 1] = (1 - np.cos(alpha)) * sy * sz + sx * np.sin(alpha)
    Q[2, 2] = np.cos(alpha) + (1 - np.cos(alpha)) * sz ** 2

def find_i_SSO(zp,za):
    R = 6378.0
    a = R + (zp + za) / 2
    e = (za - zp) / (2 * R + zp + za)
    res = root(lambda i:f(i,a,e), x0 = 97*pi/180)
    return res.x[0]*180/pi

def f(i,a,e):
    R = 6378.0 # [km]
    J2 = 1.08263e-3  # []
    mu = 3.986004418e5  # [km^3/s^2]
    O_dot = 2 * pi / (365.242 * 24 * 3600)
    return (O_dot/np.cos(i) + (3/2)*(np.sqrt(mu)*J2*R**2/((1-e**2)**2*(a)**(7/2))))