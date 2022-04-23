import numpy as np

np.set_printoptions(precision=8)


def rotation_mat(omega, I, Omega):
    P_1 = np.array([[np.cos(omega), -np.sin(omega), 0],
                    [np.sin(omega), np.cos(omega), 0],
                    [0, 0, 1]])

    P_2 = np.array([[1, 0, 0],
                    [0, np.cos(I), -np.sin(I)],
                    [0, np.sin(I), np.cos(I)]])

    P_3 = np.array([[np.cos(Omega), -np.sin(Omega), 0],
                    [np.sin(Omega), np.cos(Omega), 0],
                    [0, 0, 1]])

    P3_P2_P1 = np.dot(P_3, np.dot(P_2, P_1))

    return P3_P2_P1


def rotation_data(result_array, omega, I, Omega):
    x1, y1, z1 = result_array[1:, 0], result_array[1:, 1], result_array[1:, 2]
    x2, y2, z2 = result_array[1:, 3], result_array[1:, 4], result_array[1:, 5]

    xv1, yv1, zv1 = result_array[1:, 6], result_array[1:, 7], result_array[1:, 8]
    xv2, yv2, zv2 = result_array[1:, 9], result_array[1:, 10], result_array[1:, 11]

    m1_orb = np.array([x1, y1, z1])
    m2_orb = np.array([x2, y2, z2])
    m1_vel = np.array([xv1, yv1, zv1])
    m2_vel = np.array([xv2, yv2, zv2])

    ro1_data = np.dot(rotation_mat(omega, I, Omega), m1_orb).T
    ro2_data = np.dot(rotation_mat(omega, I, Omega), m2_orb).T
    ro1v_data = np.dot(rotation_mat(omega, I, Omega), m1_vel).T
    ro2v_data = np.dot(rotation_mat(omega, I, Omega), m2_vel).T

    ro_data = np.hstack([ro1_data, ro2_data, ro1v_data, ro2v_data])

    return ro_data



