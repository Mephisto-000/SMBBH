import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import FuncAnimation


class Plot_Result:
    def __init__(self, result_dict):
        self.no_rot_data = result_dict['no_rot_data']
        self.rot_data = result_dict['rot_data']
        self.total_energy = result_dict['total_energy']
        self.init_energy = result_dict['initial_energy']
        self.time_length = result_dict['time_length']
        self.dt = result_dict['dt']

    def plot_rk4_result(self, radius, mode="rotation_data"):
        if mode == "rotation_data":
            result_array = self.rot_data
        else:
            result_array = self.no_rot_data

        ax_len = radius / 2
        x1, y1, z1 = result_array[1:, 0], result_array[1:, 1], result_array[1:, 2]
        x2, y2, z2 = result_array[1:, 3], result_array[1:, 4], result_array[1:, 5]
        # fig = plt.figure(figsize=(30, 10))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.set_xlim3d(-ax_len, ax_len)
        ax.set_ylim3d(-ax_len, ax_len)
        ax.set_zlim3d(-ax_len, ax_len)
        ax.plot3D(x1, y1, z1, "blue", label="M1")
        ax.plot3D(x2, y2, z2, "red", label="M2")
        ax.set_xlabel('X [kpc]', fontsize=10)
        ax.set_ylabel('Y [kpc]', fontsize=10)
        ax.set_zlabel('Z [kpc]', fontsize=10)
        ax.legend(loc="upper left", fontsize=15)
        plt.show()

    def plot_total_energy(self):
        Et = self.total_energy
        plt.style.use("ggplot")
        dt_len = np.linspace(0, self.time_length, self.time_length - 1)
        plt.figure(figsize=(10, 4))
        plt.subplot()
        plt.plot(dt_len, Et, "-", color="darkblue")
        plt.title(r'$E_{total}$,' + r' ($Max(|E_{total}|)\approx$' + f'{round(max(np.abs(Et)), 4)})', fontsize=20)
        plt.xlabel("Time ($4.3*10^{15}$ years)", fontsize=10)
        plt.show()

    def plot_total_energy_divid_initE(self):
        Et = self.total_energy
        E0 = self.init_energy
        Et_E0 = Et / E0
        plt.style.use("ggplot")
        dt_len = np.linspace(0, self.time_length, self.time_length - 1)
        plt.figure(figsize=(10, 4))
        plt.subplot()
        plt.plot(dt_len, Et_E0, "-", color="darkblue")
        plt.title(r'$E_{total}/E_{0}$,' + rf' ($E_{0}\approx$ {round(E0, 4)})', fontsize=20)
        plt.xlabel("Time ($4.3*10^{15}$ years)", fontsize=10)
        plt.show()

    def plot_orbit_video(self, radius, mode="rotation_data", show_mode='plot', title=None):
        if mode == "rotation_data":
            result_array = self.rot_data[1::50, :]
        else:
            result_array = self.no_rot_data[1::50, :]
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')

        header_1 = [ax.scatter(result_array[1:, 0][0], result_array[1:, 1][0], result_array[1:, 2][0],
                               color="darkblue", marker="o", s=85, label="M1")]
        header_2 = [ax.scatter(result_array[1:, 3][0], result_array[1:, 4][0], result_array[1:, 5][0],
                               color="red", marker="o", s=85, label="M2")]

        ax.set_xlabel("X", fontsize=14)
        ax.set_ylabel("Y", fontsize=14)
        ax.set_zlabel("Z", fontsize=14)
        ax.set_title("Visualization Orbits of 2-body system\n", fontsize=16)
        ax.legend(loc="upper left", fontsize=14)

        def ani(iter, M1, M2):
            header_1[0].remove()
            header_2[0].remove()
            orb1 = ax.plot(result_array[1:, 0][:iter], result_array[1:, 1][:iter], result_array[1:, 2][:iter],
                           color="darkblue")
            orb2 = ax.plot(result_array[1:, 3][:iter], result_array[1:, 4][:iter], result_array[1:, 5][:iter],
                           color="red")
            header_1[0] = ax.scatter(result_array[1:, 0][iter], result_array[1:, 1][iter], result_array[1:, 2][iter],
                                     color="darkblue", marker="o", s=85)
            header_2[0] = ax.scatter(result_array[1:, 3][iter], result_array[1:, 4][iter], result_array[1:, 5][iter],
                                     color="red", marker="o", s=85)

            return orb1, orb2, header_1, header_2

        frames = result_array[1:, 0].shape[0]
        interval = 0.001
        anim_2b = FuncAnimation(fig, ani, frames=frames, interval=interval, repeat=False,
                                blit=False, save_count=100,
                                fargs=(header_1, header_2))

        if show_mode == 'plot':
            plt.show()
        elif show_mode == 'save':
            ani_writer = animation.writers['ffmpeg']
            writer = ani_writer(fps=10, metadata=dict(artist="Ling-Hao Lin"), bitrate=-1)
            if title is None:
                anim_2b.save(f"test.mp4", writer=writer, dpi=300)
            else:
                anim_2b.save(f"{title}.mp4", writer=writer, dpi=300)
        else:
            return anim_2b
