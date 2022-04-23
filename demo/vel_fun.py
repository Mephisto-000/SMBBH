import numpy as np

np.set_printoptions(precision=8)


def init_vel(m2, R, G):
    return np.sqrt(G*m2/R)


def init_vel_elliptical(e, init_vel_c):
    return np.sqrt((1-e)*init_vel_c**2)


def period_set(R, p=1):
    time_length = p*int(864*2*np.pi*np.sqrt(R**3)) # e = 0.5
    dt = 0.0005
    return time_length, dt


def compare_vel(v1_z, v2_z):
    vz1_vz2 = v1_z / v2_z
    return vz1_vz2

