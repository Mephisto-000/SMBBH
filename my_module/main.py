import numpy as np
from smbbh_nu import SMBBH_NU
from plot_process import Plot_Result


def potential_1(c, tmp_array, ratio_mass, r):
    ratio_term = c*(tmp_array*ratio_mass/r)
    first_term = 1/(r + r**3)
    second_term = np.arctan(r)/(r**2)
    V = ratio_term*(first_term - second_term)
    return V


if __name__ == "__main__":
    black_hole_mass = [0.5, 0.5]  # m1, m2
    angles = [np.pi/6, np.pi/4, np.pi/6]  # omega, I, Omega
    experiment_1 = SMBBH_NU(black_hole_mass, radius=0.2, eccentricity=0.5, angles=angles,
                            potential_function=potential_1)
    result_1 = experiment_1.all_result_output()

    plot1 = Plot_Result(result_dict=result_1)
    plot1.plot_rk4_result(radius=0.2)
    plot1.plot_proj_z_vel()
    plot1.plot_mv_err()
    plot1.plot_z_vel_ratio_result()
    plot1.plot_orbit_video(radius=0.2)