import numpy as np

np.set_printoptions(precision=8)


def projection_z_axis(result_array, bh_mass):
    m1, m2 = bh_mass[0], bh_mass[1]
    axis_z = np.array([0, 0, 1]).reshape((3, 1))
    m1_v = result_array[:, 6:9]
    m2_v = result_array[:, 9:12]

    ro_m1_v_zaxis = np.dot(m1_v, axis_z)
    ro_m2_v_zaxis = np.dot(m2_v, axis_z)

    v1_z = ro_m1_v_zaxis
    v2_z = ro_m2_v_zaxis

    # m1*Vz1 & m2*Vz2 : 
    m1_v1 = ro_m1_v_zaxis * m1
    m2_v2 = ro_m2_v_zaxis * m2

    return m1_v1, m2_v2, v1_z, v2_z


def mv_err(m1_v1, m2_v2):
    f = m1_v1 - m2_v2
    f = f.flatten()
    return f




