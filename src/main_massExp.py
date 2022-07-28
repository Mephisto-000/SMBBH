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

from parameters import rot_para, init_para, mass_exp



def main():
    os.system("clear")
    G = init_para["gravitational_constant"]
    R = init_para["distance"]
    e = init_para["eccentricity"]
    gal_m = init_para["galactic_mass"]


    omega = rot_para["omega"]
    I = rot_para["I"]
    Omega = rot_para["Omega"]

    m1_1 = mass_exp["m1_1"]
    m1_2 = mass_exp["m1_2"]
    m1_3 = mass_exp["m1_3"]
    m1_4 = mass_exp["m1_4"]
    m1_5 = mass_exp["m1_5"]
    m2 = mass_exp["m2"]

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

    mass_set1 = [m1_1, m2]
    mass_set2 = [m1_2, m2]
    mass_set3 = [m1_3, m2]
    mass_set4 = [m1_4, m2]
    mass_set5 = [m1_5, m2]
    mass_sets = [mass_set1, mass_set2, mass_set3, mass_set4, mass_set5]
    
    exp_set = []
    exp_rot_set = []
    exp_Vz_ratio = []
    
    for m_set in tqdm(mass_sets):
        exp = rk4(two_body_func, barycentric, init_array, total_eqs, 
                    time_length, dt, 
                    m_set, gal_m, G)
        exp_set.append(exp)
        pass
    
    for exp in exp_set:
        rot_exp = rotation_data(exp, omega, I, Omega)
        # plot_rk4_result(rot_exp, R)  # Check orbit is stable
        exp_rot_set.append(rot_exp)
    
    # plot_orbit_video(exp_rot_set[4], R)  # Output the video
    
    for rot_set, m_set in zip(exp_rot_set, mass_sets):
        m1v1, m2v2, Vz1, Vz2 = projection_z_axis(rot_set, m_set)
        Vz_ratio = compare_vel(Vz1, Vz2)
        exp_Vz_ratio.append(Vz_ratio)
    

    total_time = time_length*dt
    dt_len = np.linspace(0, total_time, time_length-1)
    plt.style.use("ggplot")
    plt.figure(figsize=(10, 5))
    plt.subplot()
    for ratio, m_set in zip(exp_Vz_ratio, mass_sets):
        plt.plot(dt_len, ratio, "-", label=f"M1 = {m_set[0]} $(M_\odot)$")
    plt.title(r"$\frac{V_{z}1}{V_{z}2}$, fix M2", fontsize=20)
    plt.xlabel("Time ($4.3\cdot 10^{15}$ years)", fontsize=10)
    
    ax = plt.gca()
    y_locator = plt.MultipleLocator(0.1)
    ax.yaxis.set_major_locator(y_locator)

    plt.legend(loc="lower right", fontsize=10)
    plt.show()


if __name__ == "__main__":
    main()
