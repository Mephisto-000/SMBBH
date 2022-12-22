import numpy as np


class SMBBH_NU:
    def __init__(self, black_holes_mass, radius, eccentricity, angles, potential_function, mass_ratio=None):
        self.__G = 1
        self.__galaxy_mass = (1e+12 / 4e+6) * 0.5
        self.period = 10
        self.__eq_amount = 6
        self.__time_length = 1000
        self.__dt = 0.0005

        self.bh_mass = black_holes_mass
        self.R = radius
        self.e = eccentricity
        self.init_c_vel = np.sqrt(self.__G*self.bh_mass[1] / self.R)
        self.init_e_vel = np.sqrt((1 - self.e)*self.init_c_vel**2)
        self.init_array = np.array([self.R, 0.0, 0.0,
                                    0.0, self.init_e_vel, 0.0])
        self.angles = angles
        self.rot_mat = np.identity(3)
        self.mass_ratio = mass_ratio
        self.poten = potential_function
        self.result_array = np.ones((self.__time_length, self.__eq_amount))
        self.result_array_barycentric_sys = np.ones((self.__time_length, self.__eq_amount*2))
        self.rotation_result = None
        self.momentum_list = []
        self.momentum_err = 0
        self.projection_z_velocity_list = []
        self.projection_z_velocity_ratio = []

    @staticmethod
    def __two_body_system(eq_index, dt_step, tmp_array, mu):
        r = np.sqrt(tmp_array[0]**2 + tmp_array[1]**2 + tmp_array[2]**2)

        if eq_index == 0:
            return tmp_array[3]
        elif eq_index == 1:
            return tmp_array[4]
        elif eq_index == 2:
            return tmp_array[5]

        elif eq_index == 3:
            return -mu * tmp_array[0] / (r ** 3)
        elif eq_index == 4:
            return -mu * tmp_array[1] / (r ** 3)
        elif eq_index == 5:
            return -mu * tmp_array[2] / (r ** 3)
        else:
            return 0

    def __barycentric(self, eq_index, dt_step, tmp_array, bh_mass, gal_mass, ratio_change):
        if ratio_change is None:
            ratio_change = [0, 0]
        m1, m2 = bh_mass[0], bh_mass[1]
        ratio_m1 = (m2 + ratio_change[0]) / (m1 + m2)
        ratio_m2 = (m1 - ratio_change[1]) / (m1 + m2)

        r1 = np.sqrt(tmp_array[0] ** 2 + tmp_array[1] ** 2 + tmp_array[2] ** 2) * ratio_m1
        r2 = np.sqrt(tmp_array[0] ** 2 + tmp_array[1] ** 2 + tmp_array[2] ** 2) * ratio_m2
        c = gal_mass * 2 / np.pi

        if eq_index == 0:
            return tmp_array[0] * ratio_m1
        elif eq_index == 1:
            return tmp_array[1] * ratio_m1
        elif eq_index == 2:
            return tmp_array[2] * ratio_m1

        elif eq_index == 3:
            return -tmp_array[0] * ratio_m2
        elif eq_index == 4:
            return -tmp_array[1] * ratio_m2
        elif eq_index == 5:
            return -tmp_array[2] * ratio_m2

        elif eq_index == 6:
            return tmp_array[3] * ratio_m1 + self.poten(c, tmp_array[0], ratio_m1, r1)
        elif eq_index == 7:
            return tmp_array[4] * ratio_m1 + self.poten(c, tmp_array[1], ratio_m1, r1)
        elif eq_index == 8:
            return tmp_array[5] * ratio_m1 + self.poten(c, tmp_array[2], ratio_m1, r1)

        elif eq_index == 9:
            return -tmp_array[3] * ratio_m2 + self.poten(c, -tmp_array[0], ratio_m2, r2)
        elif eq_index == 10:
            return -tmp_array[4] * ratio_m2 + self.poten(c, -tmp_array[1], ratio_m2, r2)
        elif eq_index == 11:
            return -tmp_array[5] * ratio_m2 + self.poten(c, -tmp_array[2], ratio_m2, r2)

    def rk4_process(self):
        m1, m2 = self.bh_mass
        mu = self.__G*(m1 + m2)

        k1, k2 = np.ones(self.__eq_amount), np.ones(self.__eq_amount)
        k3, k4 = np.ones(self.__eq_amount), np.ones(self.__eq_amount)
        v1, v2, v3 = np.ones(self.__eq_amount), np.ones(self.__eq_amount), np.ones(self.__eq_amount)
        time_tmp = 1
        dt_step = 0

        for index in range(self.__eq_amount):
            self.result_array[0, index] = self.init_array[index]

        while time_tmp < self.__time_length:
            for index in range(self.__eq_amount):
                k1[index] = self.__dt*self.__two_body_system(index, dt_step, self.result_array[time_tmp - 1], mu)
                v1[index] = self.result_array[time_tmp - 1][index] + 0.5*k1[index]
            for index in range(self.__eq_amount):
                k2[index] = self.__dt*self.__two_body_system(index, dt_step + 0.5*self.__dt, v1, mu)
                v2[index] = self.result_array[time_tmp - 1][index] + 0.5*k2[index]
            for index in range(self.__eq_amount):
                k3[index] = self.__dt*self.__two_body_system(index, dt_step + 0.5*self.__dt, v2, mu)
                v3[index] = self.result_array[time_tmp - 1][index] + k3[index]
            for index in range(self.__eq_amount):
                k4[index] = self.__dt*self.__two_body_system(index, dt_step + self.__dt, v3, mu)

            for index in range(self.__eq_amount):
                estimate_term = (k1[index] + 2*k2[index] + 2*k3[index] + k4[index]) / 6.0
                self.result_array[time_tmp][index] = self.result_array[time_tmp - 1][index] + estimate_term

            # Barycentric System :
            for index in range(self.__eq_amount*2):
                self.result_array_barycentric_sys[time_tmp][index] = self.__barycentric(index, dt_step,
                                                                                        self.result_array[time_tmp],
                                                                                        self.bh_mass,
                                                                                        self.__galaxy_mass,
                                                                                        self.mass_ratio)

            time_tmp += 1
            dt_step += self.__dt
        print("rk4 process done !")

    def rotation_mat(self):
        omega, I, Omega = self.angles
        p_1 = np.array([[np.cos(omega), -np.sin(omega), 0],
                        [np.sin(omega), np.cos(omega), 0],
                        [0, 0, 1]])
        p_2 = np.array([[1, 0, 0],
                        [0, np.cos(I), -np.sin(I)],
                        [0, np.sin(I), np.cos(I)]])
        p_3 = np.array([[np.cos(Omega), -np.sin(Omega), 0],
                        [np.sin(Omega), np.cos(Omega), 0],
                        [0, 0, 1]])

        p3p2p1 = np.dot(p_3, np.dot(p_2, p_1))
        self.rot_mat = p3p2p1

    def rotation_data(self):
        result_array = self.result_array_barycentric_sys

        x1, y1, z1 = result_array[1:, 0], result_array[1:, 1], result_array[1:, 2]
        x2, y2, z2 = result_array[1:, 3], result_array[1:, 4], result_array[1:, 5]

        xv1, yv1, zv1 = result_array[1:, 6], result_array[1:, 7], result_array[1:, 8]
        xv2, yv2, zv2 = result_array[1:, 9], result_array[1:, 10], result_array[1:, 11]

        m1_orb = np.array([x1, y1, z1])
        m2_orb = np.array([x2, y2, z2])
        m1_vel = np.array([xv1, yv1, zv1])
        m2_vel = np.array([xv2, yv2, zv2])

        self.rotation_mat()
        ro1_data = np.dot(self.rot_mat, m1_orb).T
        ro2_data = np.dot(self.rot_mat, m2_orb).T
        ro1v_data = np.dot(self.rot_mat, m1_vel).T
        ro2v_data = np.dot(self.rot_mat, m2_vel).T

        ro_data = np.hstack([ro1_data, ro2_data, ro1v_data, ro2v_data])
        self.rotation_result = ro_data

    def projection_z_axis(self):
        m1, m2 = self.bh_mass
        axis_z = np.array([0, 0, 1]).reshape((3, 1))
        m1_v = self.rotation_result[:, 6:9]
        m2_v = self.rotation_result[:, 9:12]

        ro_m1_v_zaxis = np.dot(m1_v, axis_z)
        ro_m2_v_zaxis = np.dot(m2_v, axis_z)

        v1_z = ro_m1_v_zaxis
        v2_z = ro_m2_v_zaxis

        # m1*Vz1 & m2*Vz2 :
        m1_v1 = ro_m1_v_zaxis * m1
        m2_v2 = ro_m2_v_zaxis * m2

        self.momentum_list = [m1_v1, m2_v2]
        self.projection_z_velocity_list = [v1_z, v2_z]
        self.projection_z_velocity_ratio = v1_z / v2_z

    def mv_err(self):
        f = self.momentum_list[0] - self.momentum_list[1]
        f = f.flatten()
        self.momentum_err = f

    # Output all result :
    def all_result_output(self):
        print("Begin Simulation of Supermassive Binary Black Holes : ")
        self.rk4_process()
        self.rotation_data()
        self.projection_z_axis()
        self.mv_err()
        print("All Processes Done !")

        all_results = {'no_rot_data': self.result_array_barycentric_sys,
                       'rot_data': self.rotation_result,
                       'momentum_list': self.momentum_list,
                       'momentum_err': self.momentum_err,
                       'projection_z_velocity_list': self.projection_z_velocity_list,
                       'projection_z_velocity_ratio': self.projection_z_velocity_ratio,
                       'time_length': self.__time_length,
                       'dt': self.__dt}

        return all_results


