import numpy as np


class SMBBH_NU:
    def __init__(self, black_holes_mass, time_length=5*10**3, dt=5*(10**(-4)), mass_ratio=None, **kwargs):
        self.__G = 1
        self.period = 10
        self.__eq_amount = 6
        self.time_length = time_length    # 1000
        self.dt = dt                      # 0.0005

        self.bh_mass = [m_i for m_i in black_holes_mass]
        self.mass_ratio = mass_ratio
        self.c = kwargs['constant_c']
        self.r = kwargs['radius']
        self.e = kwargs['eccentricity']
        self.angles = kwargs['angles']
        self.rot_mat = np.identity(3)
        self.potential = kwargs['potential_function']

        self.init_c_vel = np.sqrt(self.__G*self.bh_mass[1] / self.r)
        self.init_e_vel = np.sqrt((1 - self.e)*(self.init_c_vel**2))

        self.init_single_bh_array = np.array([self.r, 0.0, 0.0,
                                              0.0, self.init_e_vel, 0.0])
        self.result_array = np.ones((self.time_length, self.__eq_amount))
        self.result_array_barycentric_sys = np.ones((self.time_length, self.__eq_amount*2))

        self.rotation_result = None
        self.energy = None
        self.energy_0 = None

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

    def __barycentric(self, eq_index, dt_step, tmp_array, bh_mass, constant_c, ratio_change):
        if ratio_change is None:
            ratio_change = [0, 0]
        m1, m2 = bh_mass[0], bh_mass[1]
        ratio_m1 = (m2 + ratio_change[0]) / (m1 + m2)
        ratio_m2 = (m1 - ratio_change[1]) / (m1 + m2)

        r1 = np.sqrt(tmp_array[0] ** 2 + tmp_array[1] ** 2 + tmp_array[2] ** 2) * ratio_m1
        r2 = np.sqrt(tmp_array[0] ** 2 + tmp_array[1] ** 2 + tmp_array[2] ** 2) * ratio_m2

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
            return tmp_array[3] * ratio_m1 + self.potential(constant_c, tmp_array[0], ratio_m1, r1)
        elif eq_index == 7:
            return tmp_array[4] * ratio_m1 + self.potential(constant_c, tmp_array[1], ratio_m1, r1)
        elif eq_index == 8:
            return tmp_array[5] * ratio_m1 + self.potential(constant_c, tmp_array[2], ratio_m1, r1)

        elif eq_index == 9:
            return -tmp_array[3] * ratio_m2 + self.potential(constant_c, -tmp_array[0], ratio_m2, r2)
        elif eq_index == 10:
            return -tmp_array[4] * ratio_m2 + self.potential(constant_c, -tmp_array[1], ratio_m2, r2)
        elif eq_index == 11:
            return -tmp_array[5] * ratio_m2 + self.potential(constant_c, -tmp_array[2], ratio_m2, r2)

    def rk4_process(self):
        m1, m2 = self.bh_mass
        mu = self.__G*(m1 + m2)

        k1, k2 = np.ones(self.__eq_amount), np.ones(self.__eq_amount)
        k3, k4 = np.ones(self.__eq_amount), np.ones(self.__eq_amount)
        v1, v2, v3 = np.ones(self.__eq_amount), np.ones(self.__eq_amount), np.ones(self.__eq_amount)
        time_tmp = 1
        dt_step = 0

        for index in range(self.__eq_amount):
            self.result_array[0, index] = self.init_single_bh_array[index]

        while time_tmp < self.time_length:
            for index in range(self.__eq_amount):
                k1[index] = self.dt*self.__two_body_system(index, dt_step, self.result_array[time_tmp - 1], mu)
                v1[index] = self.result_array[time_tmp - 1][index] + 0.5*k1[index]
            for index in range(self.__eq_amount):
                k2[index] = self.dt*self.__two_body_system(index, dt_step + 0.5*self.dt, v1, mu)
                v2[index] = self.result_array[time_tmp - 1][index] + 0.5*k2[index]
            for index in range(self.__eq_amount):
                k3[index] = self.dt*self.__two_body_system(index, dt_step + 0.5*self.dt, v2, mu)
                v3[index] = self.result_array[time_tmp - 1][index] + k3[index]
            for index in range(self.__eq_amount):
                k4[index] = self.dt*self.__two_body_system(index, dt_step + self.dt, v3, mu)

            for index in range(self.__eq_amount):
                estimate_term = (k1[index] + 2*k2[index] + 2*k3[index] + k4[index]) / 6.0
                self.result_array[time_tmp][index] = self.result_array[time_tmp - 1][index] + estimate_term

            # Barycentric System :
            for index in range(self.__eq_amount*2):
                self.result_array_barycentric_sys[time_tmp][index] = self.__barycentric(index, dt_step,
                                                                                        self.result_array[time_tmp],
                                                                                        self.bh_mass,
                                                                                        self.c,
                                                                                        self.mass_ratio)

            time_tmp += 1
            dt_step += self.dt
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

        x1, y1, z1 = result_array[0:, 0], result_array[0:, 1], result_array[0:, 2]
        x2, y2, z2 = result_array[0:, 3], result_array[0:, 4], result_array[0:, 5]

        xv1, yv1, zv1 = result_array[0:, 6], result_array[0:, 7], result_array[0:, 8]
        xv2, yv2, zv2 = result_array[0:, 9], result_array[0:, 10], result_array[0:, 11]

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

    def cal_total_energy(self):
        result_array = self.rotation_result
        mass_sum = np.sum(self.bh_mass)

        x1, y1, z1 = result_array[0:, 0], result_array[0:, 1], result_array[0:, 2]
        x2, y2, z2 = result_array[0:, 3], result_array[0:, 4], result_array[0:, 5]

        r1 = np.sqrt(x1**2 + y1**2 + z1**2)
        r2 = np.sqrt(x2**2 + y2**2 + z2**2)

        xv1, yv1, zv1 = result_array[0:, 6], result_array[0:, 7], result_array[0:, 8]
        xv2, yv2, zv2 = result_array[0:, 9], result_array[0:, 10], result_array[0:, 11]

        E1 = 0.5*(xv1**2 + yv1**2 + zv1**2) - (mass_sum/r1) - self.c*(np.arctan(r1)/r1)
        E2 = 0.5*(xv2**2 + yv2**2 + zv2**2) - (mass_sum/r2) - self.c*(np.arctan(r2)/r2)
        total_E = E1[1:] + E2[1:]

        self.energy = total_E
        self.energy_0 = total_E[0]

    # Output all result :
    def run(self):
        print("Begin Simulation of Supermassive Binary Black Holes : ")
        self.rk4_process()
        self.rotation_data()
        self.cal_total_energy()
        print("All Processes Done !")

        all_results = {'no_rot_data': self.result_array_barycentric_sys,
                       'rot_data': self.rotation_result,
                       'total_energy': self.energy,
                       'initial_energy': self.energy_0,
                       'time_length': self.time_length,
                       'dt': self.dt}

        return all_results
