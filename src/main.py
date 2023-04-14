import os
import numpy as np
from smbbh_nu import SMBBH_NU
from plot_process import Plot_Result


def potential_1(constant_c, comp_vector, r):
    c_term = (constant_c*comp_vector) / r
    first_term = 1/(r + r**3)
    second_term = np.arctan(r)/(r**2)
    V = c_term*(first_term - second_term)
    return V


if __name__ == "__main__":
    black_hole_mass = [0.5, 0.5]  # m1(0.5), m2(0.5)
    t_0 = 0
    t_f = 50
    radius = 2  # 2
    eccentricity = 0.7  # 0.7
    angles = [np.pi/6, np.pi/4, np.pi/6]  # omega, I, Omega
    # angles = [0.0, 0.0, 0.0]  # omega, I, Omega
    constant_c = 0.9  # 1.0, 0.99, 0.9, 0.8

    experiment_1 = SMBBH_NU(black_hole_mass,
                            t0=t_0,
                            tf=t_f,
                            constant_c=constant_c,
                            radius=radius,
                            eccentricity=eccentricity,
                            angles=angles,
                            potential_function=potential_1)

    exp_dict = experiment_1.run()

    while True:
        os.system("clear")
        print()
        print("To display the trajectory chart for the supermassive binary black hole orbit, please press 'c1'")
        print("To display the total energy change chart for the supermassive binary black holes, please press 'c2'")
        print("To display the total energy change divided by the initial total energy chart for the "
              "supermassive binary black holes, please press 'c3'")
        print("To display the video of the trajectory chart for the supermassive binary black hole "
              "orbit, please press 'c4'")
        print("To save the video of the trajectory chart for the supermassive binary black hole "
              "orbit, please press 'c5'")
        print("To exit, please press 'q'")
        print()
        print()

        test_command = str(input("Please input the test case: "))
        if test_command == "q":
            break

        rot_command = str(input("Would you like to rotate the binary black hole orbit trajectory in 3D? [y]/n : "))
        rot_mod = "rotation" if rot_command == "y" else "no_rotation"
        plot = Plot_Result(exp_dict, radius, t_0, t_f)

        if test_command == "c1":
            plot.plot_rk4_result(mode=rot_mod)
        elif test_command == "c2":
            plot.plot_total_energy(mode=rot_mod)
        elif test_command == "c3":
            plot.plot_total_energy_divid_initE(mode=rot_mod)
        elif test_command == "c4":
            plot.plot_orbit_video(mode=rot_mod)
        elif test_command == "c5":
            plot.plot_orbit_video(mode=rot_mod, show_mode="save")
        else:
            print("Input Error! Please input again.")
