import numpy as np



def cent_pos(data, mass_list):
    m1, m2 = mass_list[0], mass_list[1]
    x1, y1, z1 = data[1:, 0], data[1:, 1], data[1:, 2]
    x2, y2, z2 = data[1:, 3], data[1:, 4], data[1:, 5]

    mx1, my1, mz1 = m1*x1, m1*y1, m1*z1
    mx2, my2, mz2 = m2*x2, m2*y2, m2*z2

    c_m_x = (mx1 + mx2) / (m1 + m2)
    c_m_y = (my1 + my2) / (m1 + m2)
    c_m_z = (mz1 + mz2) / (m1 + m2)

    return c_m_x, c_m_y, c_m_z





