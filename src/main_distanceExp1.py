import os
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from vel_fun import init_vel, init_vel_elliptical, period_set, compare_vel
from rk4 import rk4
from eq_sys import two_body_func, barycentric
from rotation_fun import rotation_data
from proj_Z import projection_z_axis, mv_err
from plot_fun import plot_rk4_result, plot_proj_z_result, \
                     plot_mv_err_result, plot_proj_z_vel, plot_vel_compare_result, \
                     plot_orbit_video

from parameters import rot_para, init_para, distance_exp1


def main():
    os.system("clear")
    G = init_para["gravitational_constant"]
    R = init_para["distance"]
    e = init_para["eccentricity"]
    gal_m = init_para["galactic_mass"]
    m1 = init_para["black_hole1_mass"]
    m2 = init_para["black_hole2_mass"]
    m_set = [m1, m2]

    omega = rot_para["omega"]
    I = rot_para["I"]
    Omega = rot_para["Omega"]

    r1_1 = distance_exp1["r_1"]
    r1_2 = distance_exp1["r_2"]
    r1_3 = distance_exp1["r_3"]
    r1_4 = distance_exp1["r_4"]
    r1_5 = distance_exp1["r_5"]
    r1_exp_set = [r1_1, r1_2, r1_3, r1_4, r1_5]

    # initial value : 
    init_v_circle = init_vel(m2, R, G)
    init_v = init_vel_elliptical(e, init_v_circle)
    init_array = np.array([R, 0.0, 0.0, 
                           0.0, init_v, 0.0])

    # rk4 equations (orbit xyz and velocity xyz) using : 
    total_eqs = 6

    # Time setting : 
    time_length = 2000
    dt = 0.0005
    print(np.format_float_scientific(time_length*dt*4.3e+15), " (years)")

    exp_set = []
    exp_rot_set = []
    exp_Vz_ratio = []

    for r1_set in tqdm(r1_exp_set):
        exp = rk4(two_body_func, barycentric, init_array, total_eqs, 
                    time_length, dt, 
                    m_set, gal_m, G, r1_set)
        exp_set.append(exp)

    for exp in exp_set:
        rot_exp = rotation_data(exp, omega, I, Omega)
        # plot_rk4_result(rot_exp, R)  # Check orbit is stable
        exp_rot_set.append(rot_exp)


    for rot_set in exp_rot_set:
        m1v1, m2v2, Vz1, Vz2 = projection_z_axis(rot_set, m_set)
        Vz_ratio = compare_vel(Vz1, Vz2)
        exp_Vz_ratio.append(Vz_ratio)

    

    total_time = time_length*dt
    dt_len = np.linspace(0, total_time, time_length-1)
    plt.style.use("ggplot")
    plt.figure(figsize=(10, 5))
    plt.subplot()
    for ratio, r1_set in zip(exp_Vz_ratio, r1_exp_set):
        plt.plot(dt_len, ratio, "-", 
                label=f"$\epsilon$ = {round(r1_set[0]*0.1, 4)} (kpc)")
    plt.title(r"$\frac{V_{z}1}{V_{z}2}$, $r_{1}+\epsilon$", fontsize=20)
    plt.xlabel("Time ($4.3\cdot 10^{15}$ years)", fontsize=10)
    
    ax = plt.gca()
    y_locator = plt.MultipleLocator(0.1)
    ax.yaxis.set_major_locator(y_locator)

    plt.legend(loc="lower right", fontsize=10)
    plt.show()


if __name__ == "__main__":
    main()
    