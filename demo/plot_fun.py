import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import FuncAnimation


def plot_rk4_result(result_array, R):
    ax_len = R / 2
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
    # plt.close()


def plot_proj_z_result(v1_m1, v2_m2, time_length, dt):
    plt.style.use("ggplot")
    total_time = time_length*dt
    dt_len = np.linspace(0, total_time, time_length-1)
    plt.plot(dt_len, v1_m1, "-", color="darkblue", label=r"$m_{1}\cdot V_{z}1$")
    plt.plot(dt_len, v2_m2, "-", color="red", label=r"$m_{2}\cdot V_{z}2$")
    plt.legend(fontsize=20)
    plt.xlabel("Time ($4.3*10^{15}$ years)", fontsize=20)
    plt.legend(loc="lower right", fontsize=20)
    plt.show()
    # plt.close()


def plot_proj_z_vel(v1_z, v2_z, time_length, dt):
    plt.style.use("ggplot")
    total_time = time_length*dt
    dt_len = np.linspace(0, total_time, time_length-1)
    plt.figure()
    plt.subplot()
    plt.plot(dt_len, v1_z, "-", color="darkblue", label=r"$V_{z}1$")
    plt.plot(dt_len, v2_z, "-", color="red", label=r"$V_{z}2$")
    plt.legend(fontsize=15)
    plt.xlabel("Time ($4.3*10^{15}$ years)", fontsize=10)
    plt.legend(loc="lower right", fontsize=10)
    plt.show()
    # plt.close()


def plot_mv_err_result(f, time_length, dt):
    total_time = time_length*dt
    dt_len = np.linspace(0, total_time, len(f))
    plt.style.use("ggplot")
    plt.plot(dt_len, f, "-", color="red", label="$m_{1}*V_{z}1 - m_{2}*V_{z}2$")
    plt.axhline(y=0, color="green", linestyle="-.", label="$y=0$")
    plt.xlabel("Time ($4.3\cdot 10^{15}$ years)", fontsize=20)
    plt.legend(loc="lower right", fontsize=20)
    plt.show()
    # plt.close()


def plot_vel_compare_result(vz1_vz2, time_length, dt, title=None):
    total_time = time_length*dt
    dt_len = np.linspace(0, total_time, len(vz1_vz2))
    plt.style.use("ggplot")
    plt.figure(figsize=(10, 5))
    plt.subplot()
    plt.plot(dt_len, vz1_vz2, "-", color="red")
    if title:
        if title[0] == "m":
            plt.title(r"$\frac{V_{z}1}{V_{z}2}, $" + f"M1={title[1]}, M2=0.5", fontsize=20)
        elif title[0] == "r":
            plt.title(r"$\frac{V_{z}1}{V_{z}2}, $" + f"$r_1$ = (1/2)(1+{title[1][0]})$\cdot$0.2, "
                                                     f"$r_2$ = (1/2)(1-{title[1][1]})$\cdot$0.2"
                      , fontsize=20)
        else:
            print("Error !")
    else:
        plt.title(r"$\frac{V_{z}1}{V_{z}2}$", fontsize=20)
    plt.xlabel("Time ($4.3\cdot 10^{15}$ years)", fontsize=10)
    plt.show()
    # plt.close()


def plot_orbit_video(result_array, R):
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
        orb2 = ax.plot(result_array[1:, 3][:iter], result_array[1:, 4][:iter], result_array[1:, 5][:iter], color="red")
        header_1[0] = ax.scatter(result_array[1:, 0][iter], result_array[1:, 1][iter], result_array[1:, 2][iter],
                                 color="darkblue", marker="o", s=85)
        header_2[0] = ax.scatter(result_array[1:, 3][iter], result_array[1:, 4][iter], result_array[1:, 5][iter],
                                 color="red", marker="o", s=85)

        return orb1, orb2, header_1, header_2

    anim_2b = FuncAnimation(fig, ani, frames=result_array[1:, 0].shape[0], interval=0.2, repeat=False, blit=False,
                            fargs=(header_1, header_2))

    ani_writer = animation.writers['ffmpeg']
    writer = ani_writer(fps=100, metadata=dict(artist="Ling-Hao Lin"), bitrate=-1)
    anim_2b.save(f"test.mp4", writer=writer, dpi=300)


