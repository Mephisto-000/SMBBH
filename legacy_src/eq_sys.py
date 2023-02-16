import numpy as np

np.set_printoptions(precision=8)


def two_body_func(eq_index, dt_step, tmp_array, mu):
    r = np.sqrt(tmp_array[0]**2 + tmp_array[1]**2 + tmp_array[2]**2)

    if eq_index == 0:
        return tmp_array[3]
    elif eq_index == 1:
        return tmp_array[4]
    elif eq_index == 2:
        return tmp_array[5]

    elif eq_index == 3:
        return -mu*tmp_array[0]/(r**3)
    elif eq_index == 4:
        return -mu*tmp_array[1]/(r**3)
    elif eq_index == 5:
        return -mu*tmp_array[2]/(r**3)

    else:
        return 0


def poten(c, tmp_array, ratio_mass, r):
    ratio_term = c*(tmp_array*ratio_mass/r)
    first_term = 1/(r + r**3)
    second_term = np.arctan(r)/(r**2)

    V = ratio_term*(first_term - second_term)
    return V


def barycentric(eq_index, dt_step, tmp_array, bh_mass, gal_mass, ratio_change=None):
    if ratio_change is None:
        ratio_change = [0, 0]
    m1, m2 = bh_mass[0], bh_mass[1]
    ratio_m1 = (m2 + ratio_change[0])/(m1 + m2)
    ratio_m2 = (m1 - ratio_change[1])/(m1 + m2)

    r1 = np.sqrt(tmp_array[0]**2 + tmp_array[1]**2 + tmp_array[2]**2)*ratio_m1
    r2 = np.sqrt(tmp_array[0]**2 + tmp_array[1]**2 + tmp_array[2]**2)*ratio_m2
    c = gal_mass*2/np.pi

    if eq_index == 0:
        return tmp_array[0]*ratio_m1
    elif eq_index == 1:
        return tmp_array[1]*ratio_m1
    elif eq_index == 2:
        return tmp_array[2]*ratio_m1

    elif eq_index == 3:
        return -tmp_array[0]*ratio_m2
    elif eq_index == 4:
        return -tmp_array[1]*ratio_m2
    elif eq_index == 5:
        return -tmp_array[2]*ratio_m2

    elif eq_index == 6:
        return tmp_array[3]*ratio_m1 + poten(c, tmp_array[0], ratio_m1, r1)
    elif eq_index == 7:
        return tmp_array[4]*ratio_m1 + poten(c, tmp_array[1], ratio_m1, r1)
    elif eq_index == 8:
        return tmp_array[5]*ratio_m1 + poten(c, tmp_array[2], ratio_m1, r1)

    elif eq_index == 9:
        return -tmp_array[3]*ratio_m2 + poten(c, -tmp_array[0], ratio_m2, r2)
    elif eq_index == 10:
        return -tmp_array[4]*ratio_m2 + poten(c, -tmp_array[1], ratio_m2, r2)
    elif eq_index == 11:
        return -tmp_array[5]*ratio_m2 + poten(c, -tmp_array[2], ratio_m2, r2)





