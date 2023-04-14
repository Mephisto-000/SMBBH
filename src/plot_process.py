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

    def plot_rk4_result(self, mode="rotation"):
        if mode == "rotation":
            data_dict = self.rot_data
            r1 = data_dict['p1_orbit']
            r2 = data_dict['p2_orbit']

        elif mode == "no_rotation":
            data_dict = self.no_rot_data
            r1 = data_dict['p1_orbit']
            r2 = data_dict['p2_orbit']
        else:
            print("Test mode")
            data_dict = self.no_rot_data
            r1 = data_dict['p1_orbit']
            r2 = data_dict['p2_orbit']

        ax_len = self.radius * 0.75
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.set_xlim3d(-ax_len, ax_len)
        ax.set_ylim3d(-ax_len, ax_len)
        ax.set_zlim3d(-ax_len, ax_len)
        ax.plot(r1[:, 0], r1[:, 1], r1[:, 2], color='darkblue', alpha=0.4,  label="p1 orbit")
        ax.plot(r2[:, 0], r2[:, 1], r2[:, 2], color='red', alpha=0.4, label="p2 orbit")
        ax.plot(r1[0, 0], r1[0, 1], r1[0, 2], "bo", label="p1 orbit (t=0)")
        ax.plot(r2[0, 0], r2[0, 1], r2[0, 2], "ro", label="p2 orbit (t=0)")
        ax.plot(0, 0, 0, 'ko', label="COG")
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
        plt.figure(figsize=(10, 4))
        plt.subplot()
        plt.plot(self.time_length, total_energy, "-", color="darkblue")
        plt.title(r'$E_{total}$,' + r' ($Max(|E_{total}|)\approx$' + f'{round(max(np.abs(total_energy)), 4)})', fontsize=20)
        plt.xlabel("Time ($4.3*10^{15}$ years)", fontsize=10)
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

        Et_divid_E0 = total_energy/initial_energy
        plt.style.use("ggplot")
        plt.figure(figsize=(10, 4))
        plt.subplot()
        plt.plot(self.time_length, Et_divid_E0, "-", color="darkblue")
        plt.title(r'$E_{total}/E_{0}$,' + rf' ($E_{0}\approx$ {round(initial_energy, 4)})', fontsize=20)
        plt.xlabel("Time ($4.3*10^{15}$ years)", fontsize=10)
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
        ax_len = self.radius * 0.75
        ax.set_xlim3d(-ax_len, ax_len)
        ax.set_ylim3d(-ax_len, ax_len)
        ax.set_zlim3d(-ax_len, ax_len)

        header_1 = [ax.scatter(r1[:, 0][0], r1[:, 1][0], r1[:, 2][0],
                               color="darkblue", marker="o", s=85, label="M1")]
        header_2 = [ax.scatter(r2[:, 0][0], r2[:, 1][0], r2[:, 2][0],
                               color="red", marker="o", s=85, label="M2")]

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

