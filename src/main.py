import os
import numpy as np
from smbbh_nu import SMBBH_NU
from plot_process import Plot_Result


def potential_1(constant_c, tmp_array, ratio_mass, r):
    ratio_term = constant_c*(tmp_array*ratio_mass/r)
    first_term = 1/(r + r**3)
    second_term = np.arctan(r)/(r**2)
    V = ratio_term*(first_term - second_term)
    return V


if __name__ == "__main__":
    black_hole_mass = [0.5, 0.5]  # m1(0.5), m2(0.5)
    radius = 0.5
    eccentricity = 0
    angles = [np.pi/6, np.pi/4, np.pi/6]  # omega, I, Omega
    constant_c = 0  # 1000
    experiment_1 = SMBBH_NU(black_hole_mass,
                            constant_c=constant_c,
                            radius=radius,
                            eccentricity=eccentricity,
                            angles=angles,
                            potential_function=potential_1)
    result_1 = experiment_1.run()

    while True:
        test_command = str(input("Please input the test case: "))
        plot1 = Plot_Result(result_dict=result_1)
        if test_command == "c1":
            plot1.plot_rk4_result(radius=radius)
            break
        elif test_command == "c2":
            plot1.plot_total_energy()
            break
        elif test_command == "c3":
            plot1.plot_total_energy_divid_initE()
            break
        elif test_command == "c4":
            show_mode = input("Please choose show mode [plot/save] : ")
            plot1.plot_orbit_video(radius=0.05, show_mode=show_mode)
            break
        else:
            os.system("clear")
            print("Input Error! Please input again.")
