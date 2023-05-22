import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import FuncAnimation


class Plot_Result:
    def __init__(self, all_result, radius, t0, tf, tl=1000):
        self.no_rot_data = all_result['no_rot_data']
        self.rot_data = all_result['rot_data']

        self.no_rot_total_energy = all_result['no_rot_energy']
        self.no_rot_init_energy = all_result['no_rot_energy0']
        self.total_energy = all_result['rot_energy']
        self.init_energy = all_result['rot_energy0']

        self.radius = radius
        self.t0 = t0
        self.tf = tf
        self.time_length = np.linspace(t0, tf, tl)

    def __mode_choose(self, mode="rotation"):
        if mode == "rotation":
            data_dict = self.rot_data
            r1 = data_dict['p1_orbit']
            r2 = data_dict['p2_orbit']
            return r1, r2
        elif mode == "no_rotation":
            data_dict = self.no_rot_data
            r1 = data_dict['p1_orbit']
            r2 = data_dict['p2_orbit']
            return r1, r2
        else:
            print("Test mode")
            data_dict = self.no_rot_data
            r1 = data_dict['p1_orbit']
            r2 = data_dict['p2_orbit']
            return r1, r2

    def plot_rk4_result(self, mode="rotation"):
        r1, r2 = self.__mode_choose(mode)
        ax_len = self.radius
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.set_xlim3d(-ax_len, ax_len)
        ax.set_ylim3d(-ax_len, ax_len)
        ax.set_zlim3d(-ax_len, ax_len)
        ax.plot(r1[:, 0], r1[:, 1], r1[:, 2], color='darkblue', alpha=0.4,  label="P1 orbit")
        ax.plot(r2[:, 0], r2[:, 1], r2[:, 2], color='red', alpha=0.4, label="P2 orbit")
        ax.plot(r1[0, 0], r1[0, 1], r1[0, 2], "bo", label="P1 orbit (t=0)")
        ax.plot(r2[0, 0], r2[0, 1], r2[0, 2], "ro", label="P2 orbit (t=0)")
        ax.plot(0, 0, 0, 'ko', label=r"O, $(0, 0, 0)$")
        ax.legend()
        plt.show()

    def plot_total_energy(self, mode="rotation"):
        if mode == "rotation":
            total_energy = self.total_energy
        elif mode == "no_rotation":
            total_energy = self.no_rot_total_energy
        else:
            print("Test mode")
            total_energy = self.no_rot_total_energy

        plt.style.use("ggplot")
        plt.figure(figsize=(10, 5))
        plt.subplot()
        plt.plot(self.time_length, total_energy, "-", color="darkblue")
        plt.xlabel("t", fontsize=20)
        plt.ylabel(r"$E(t)$        ", fontsize=20, rotation=0)
        plt.show()

    def plot_total_energy_divid_initE(self, mode="rotation"):
        if mode == "rotation":
            total_energy = self.total_energy
            initial_energy = self.init_energy
        elif mode == "no_rotation":
            total_energy = self.no_rot_total_energy
            initial_energy = self.no_rot_init_energy
        else:
            print("Test mode")
            total_energy = self.no_rot_total_energy
            initial_energy = self.no_rot_init_energy

        if initial_energy != 0:
            Et_divid_E0 = total_energy/initial_energy
            plt.style.use("ggplot")
            plt.figure(figsize=(10, 5))
            plt.subplot()
            plt.plot(self.time_length, Et_divid_E0, "-", color="darkblue")
            plt.xlabel("t", fontsize=20)
            plt.ylabel(r"$\frac{E(t)}{E_{0}}$        ", fontsize=20, rotation=0)
            plt.show()
        else:
            print("\033[31m initial energy is 0. \033[0m")
            enter_key = input("")

    def plot_2D_plane_orbit(self, mode="rotation", plane_choose="x-y"):
        r1, r2 = self.__mode_choose(mode)
        plt.style.use("ggplot")
        plt.figure()
        plt.subplot()
        if plane_choose == "x-y":
            plt.plot(r1[:, 0], r1[:, 1], color="darkblue", label=r"$P_{1}$ orbit", alpha=0.4)
            plt.plot(r2[:, 0], r2[:, 1], color="red", label=r"$P_{2}$ orbit", alpha=0.4)
            plt.plot(0, 0, 'ko', label=r"$(0, 0)$")
            plt.xlabel(r"$\hat{X}$", fontsize=18)
            plt.ylabel(r"$\hat{Y}$", fontsize=18, rotation=0)
        elif plane_choose == "x-z":
            plt.plot(r1[:, 0], r1[:, 2], color="darkblue", label=r"$P_{1}$ orbit", alpha=0.4)
            plt.plot(r2[:, 0], r2[:, 2], color="red", label=r"$P_{2}$ orbit", alpha=0.4)
            plt.plot(0, 0, 'ko', label=r"$(0, 0)$")
            plt.xlabel(r"$\hat{X}$", fontsize=18)
            plt.ylabel(r"$\hat{Z}$", fontsize=18, rotation=0)
        elif plane_choose == "y-z":
            plt.plot(r1[:, 1], r1[:, 2], color="darkblue", label=r"$P_{1}$ orbit", alpha=0.4)
            plt.plot(r2[:, 1], r2[:, 2], color="red", label=r"$P_{2}$ orbit", alpha=0.4)
            plt.plot(0, 0, 'ko', label=r"$(0, 0)$")
            plt.xlabel(r"$\hat{Y}$", fontsize=18)
            plt.ylabel(r"$\hat{Z}$", fontsize=18, rotation=0)
        else:
            print("Error")
        plt.axis("square")
        plt.legend(bbox_to_anchor=(1.35, 1))
        plt.show()

    def plot_2D_time_orbit(self, mode="rotation", plane_choose="x"):
        r1, r2 = self.__mode_choose(mode)
        plt.style.use("ggplot")
        plt.figure(figsize=(10, 5))
        plt.subplot()
        if plane_choose == "x":
            plt.plot(self.time_length, r1[:, 0], color="darkblue", label="P1")
            plt.plot(self.time_length, r2[:, 0], color="red", label="P2")
            plt.xlabel("t", fontsize=20)
            plt.ylabel("X", fontsize=20, rotation=0)
            plt.legend(bbox_to_anchor=(1.1, 1))
            plt.show()
        elif plane_choose == "y":
            plt.plot(self.time_length, r1[:, 1], color="darkblue", label="P1")
            plt.plot(self.time_length, r2[:, 1], color="red", label="P2")
            plt.xlabel("t", fontsize=20)
            plt.ylabel("Y", fontsize=20, rotation=0)
            plt.legend(bbox_to_anchor=(1.1, 1))
            plt.show()
        elif plane_choose == "z":
            plt.plot(self.time_length, r1[:, 2], color="darkblue", label="P1")
            plt.plot(self.time_length, r2[:, 2], color="red", label="P2")
            plt.xlabel("t", fontsize=20)
            plt.ylabel("Z", fontsize=20, rotation=0)
            plt.legend(bbox_to_anchor=(1.1, 1))
            plt.show()
        else:
            print("Error")

    def plot_r_length_per_time(self, mode="rotation", particle="p1"):
        r1, r2 = self.__mode_choose(mode)
        plt.style.use("ggplot")
        plt.figure(figsize=(10, 5))
        plt.subplot()
        if particle == "p1":
            x_orb, y_orb, z_orb = r1[:, 0], r1[:, 1], r1[:, 2]
            r = np.sqrt(x_orb**2 + y_orb**2 + z_orb**2)
            plt.plot(self.time_length, r, color="darkblue")
            plt.xlabel("t", fontsize=20)
            plt.ylabel(r"$r_{1}$        ", fontsize=20, rotation=0)
        elif particle == "p2":
            x_orb, y_orb, z_orb = r2[:, 0], r2[:, 1], r2[:, 2]
            r = np.sqrt(x_orb**2 + y_orb**2 + z_orb**2)
            plt.plot(self.time_length, r, color="red")
            plt.xlabel("t", fontsize=20)
            plt.ylabel(r"$r_{2}$        ", fontsize=20, rotation=0)
        else:
            print("Error")
        plt.show()

    def plot_r12_length_per_time(self, mode="rotation"):
        r1, r2 = self.__mode_choose(mode)
        plt.style.use("ggplot")
        plt.figure(figsize=(10, 5))
        plt.subplot()
        r12_x = r2[:, 0] - r1[:, 0]
        r12_y = r2[:, 1] - r1[:, 1]
        r12_z = r2[:, 2] - r1[:, 2]
        r12_len = np.sqrt(r12_x**2 + r12_y**2 + r12_z**2)
        print(f"\n minimum length : {np.round(np.min(r12_len), 2)}")
        plt.plot(self.time_length, r12_len, color="green")
        plt.xlabel("t", fontsize=20)
        plt.ylabel(r"$r_{12}$        ", fontsize=20, rotation=0)
        plt.show()

    def plot_orbit_video(self, mode="rotation", show_mode='plot', title=None):
        if mode == "rotation":
            data_dict = self.rot_data
            r1 = data_dict['p1_orbit'][::5, :]
            r2 = data_dict['p2_orbit'][::5, :]
        elif mode == "no_rotation":
            data_dict = self.no_rot_data
            r1 = data_dict['p1_orbit'][::5, :]
            r2 = data_dict['p2_orbit'][::5, :]
        else:
            print("Test mode")
            data_dict = self.no_rot_data
            r1 = data_dict['p1_orbit']
            r2 = data_dict['p2_orbit']

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax_len = self.radius
        ax.set_xlim3d(-ax_len, ax_len)
        ax.set_ylim3d(-ax_len, ax_len)
        ax.set_zlim3d(-ax_len, ax_len)

        header_1 = [ax.scatter(r1[:, 0][0], r1[:, 1][0], r1[:, 2][0],
                               color="darkblue", marker="o", s=85, label="P1")]
        header_2 = [ax.scatter(r2[:, 0][0], r2[:, 1][0], r2[:, 2][0],
                               color="red", marker="o", s=85, label="P2")]

        ax.set_xlabel("X", fontsize=14)
        ax.set_ylabel("Y", fontsize=14)
        ax.set_zlabel("Z", fontsize=14)
        ax.set_title("Animation of the orbit trajectory for the SMBBH\n", fontsize=16)
        ax.legend(loc="upper left", fontsize=14)

        def ani(iter, M1, M2):
            header_1[0].remove()
            header_2[0].remove()
            orb1 = ax.plot(r1[:, 0][:iter], r1[:, 1][:iter], r1[:, 2][:iter],
                           color="darkblue", alpha=0.3)
            orb2 = ax.plot(r2[:, 0][:iter], r2[:, 1][:iter], r2[:, 2][:iter],
                           color="red", alpha=0.3)
            header_1[0] = ax.scatter(r1[:, 0][iter], r1[:, 1][iter], r1[:, 2][iter],
                                     color="darkblue", marker="o", s=85)
            header_2[0] = ax.scatter(r2[:, 0][iter], r2[:, 1][iter], r2[:, 2][iter],
                                     color="red", marker="o", s=85)

            return orb1, orb2, header_1, header_2

        frames = r1[:, 0].shape[0]
        interval = 0.001
        anim_2b = FuncAnimation(fig, ani, frames=frames, interval=interval, repeat=False,
                                blit=False, save_count=100,
                                fargs=(header_1, header_2))

        if show_mode == 'plot':
            plt.show()
        elif show_mode == 'save':
            ani_writer = animation.writers['ffmpeg']
            writer = ani_writer(fps=10, metadata=dict(artist="Lin-Hao Lin"), bitrate=-1)
            if title is None:
                anim_2b.save(f"test.mp4", writer=writer, dpi=300)
            else:
                anim_2b.save(f"{title}.mp4", writer=writer, dpi=300)
        else:
            return anim_2b
