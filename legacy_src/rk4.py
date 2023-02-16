import numpy as np

np.set_printoptions(precision=8)


def rk4(eq, eq_bc, init_val, eq_amo, t_len, dt, bh_mass, gal_mass, G, ratio_change=None):
    """

    :param eq: OD equation
    :param init_val: initial values
    :param eq_amo: equation amount
    :param t_len: time length
    :param dt: time step
    :return: [m1_orb_x, m1_orb_y, m1_orb_z, m2_orb_x, m2_orb_y, m2_orb_z,
              m1_vel_x, m1_vel_y, m1_vel_z, m2_vel_x, m2_vel_y, m2_vel_z]
    """
    m1, m2 = bh_mass[0], bh_mass[1]
    mu = G*(m1 + m2)
    result_array = np.ones((t_len, eq_amo))

    # Barycentric System array :
    result_array_barycentric = np.ones((t_len, eq_amo*2))

    k1, k2, k3, k4 = np.ones(eq_amo), np.ones(eq_amo), np.ones(eq_amo), np.ones(eq_amo)
    v1, v2, v3 = np.ones(eq_amo), np.ones(eq_amo), np.ones(eq_amo)
    time_tmp = 1
    dt_step = 0

    for index in range(eq_amo):
        result_array[0, index] = init_val[index]

    while time_tmp < t_len:
        for index in range(eq_amo):
            k1[index] = dt*eq(index, dt_step, result_array[time_tmp-1], mu)
            v1[index] = result_array[time_tmp-1][index] + 0.5*k1[index]
        for index in range(eq_amo):
            k2[index] = dt*eq(index, dt_step + 0.5*dt, v1, mu)
            v2[index] = result_array[time_tmp-1][index] + 0.5*k2[index]
        for index in range(eq_amo):
            k3[index] = dt*eq(index, dt_step + 0.5*dt, v2, mu)
            v3[index] = result_array[time_tmp-1][index] + k3[index]
        for index in range(eq_amo):
            k4[index] = dt*eq(index, dt_step + dt, v3, mu)

        for index in range(eq_amo):
            estimate_term = (k1[index] + 2*k2[index] + 2*k3[index] + k4[index]) / 6.0
            result_array[time_tmp][index] = result_array[time_tmp-1][index] + estimate_term

        # Barycentric System :
        for index in range(eq_amo*2):
            result_array_barycentric[time_tmp][index] = eq_bc(index, dt_step, result_array[time_tmp],
                                                              bh_mass,
                                                              gal_mass,
                                                              ratio_change)

        time_tmp += 1
        dt_step += dt
    print("Done !")

    return result_array_barycentric


