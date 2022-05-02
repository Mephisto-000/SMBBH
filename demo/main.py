import numpy as np

from vel_fun import init_vel, init_vel_elliptical, period_set, compare_vel
from rk4 import rk4
from eq_sys import two_body_func, barycentric
from rotation_fun import rotation_data
from proj_Z import projection_z_axis, mv_err
from plot_fun import plot_rk4_result, plot_proj_z_result, \
                     plot_mv_err_result, plot_proj_z_vel, plot_vel_compare_result, \
                     plot_orbit_video


"""
Mass unit   : 1.2e+12 (Sun mass)
Length unit : kpc
Time unit   : 4.3e+15 (Years)
G = 1

Mass of Galaxy (Milky Way) : ((0.8~1.5) * 1e+12).       Choose 1e+12 (Sun Mass)
Black Hole (Milky Way)     : ((4.154 +/-0.014)*1e+6).   Choose 4e+6  (Sun Mass)
"""


if __name__ == "__main__":
    m1 = 0.5  # 0.9
    m2 = 0.5  # fix
    bh_mass = [m1, m2]
    gal_mass = (1e+12 / 4e+6) * 0.5

    G = 1
    R = 0.2
    e = 0.8  # 0.8

    omega = np.pi / 6
    I = np.pi / 4
    Omega = np.pi / 6

    period = 10
    total_amount = 6

    init_v = init_vel(m2, R, G)
    init_v = init_vel_elliptical(e, init_v)
    time_length = 1000
    dt = 0.0005

    print(time_length*dt)
    init_array = np.array([R, 0.0, 0.0,
                           0.0, init_v, 0.0])

    rk4_result = rk4(two_body_func, barycentric,
                     init_array, total_amount, time_length, dt,
                     bh_mass, gal_mass, G)
    plot_rk4_result(rk4_result, R)

    rot_data = rotation_data(rk4_result, omega, I, Omega)
    plot_rk4_result(rot_data, R)


    m1_v1, m2_v2, v1_z, v2_z = projection_z_axis(rot_data, bh_mass)
    plot_proj_z_vel(v1_z, v2_z, time_length, dt)
    plot_proj_z_result(m1_v1, m2_v2, time_length, dt)

    f = mv_err(m1_v1, m2_v2)
    # plot_mv_err_result(f, time_length, dt)

    vz1_vz2 = compare_vel(v1_z, v2_z)
    plot_vel_compare_result(vz1_vz2, time_length, dt)

    # plot_orbit_video(rot_data, R)
