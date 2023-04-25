import numpy as np


class SMBBH_NU:
    def __init__(self, black_holes, t0, tf, tl=1000, **kwargs):
        self.__G = 1
        self.__eq_amount = 12
        self.__tl = tl
        self.t0 = t0
        self.tf = tf
        self.__dt = (tf - t0) / (tl - 1)
        self.__t_p = np.linspace(t0, tf, tl)

        self.bh_mass = {'m1': black_holes[0],
                        'm2': black_holes[1]}
        self.c = kwargs['constant_c']
        self.r = kwargs['radius']
        self.e = kwargs['eccentricity']
        self.angles = kwargs['angles']
        self.rot_mat = np.identity(3)
        self.potential = kwargs['potential_function']

        self.init_c_vel = np.sqrt(self.__G*self.bh_mass['m1'] / self.r)
        self.init_e_vel = np.sqrt((1 - self.e)*(self.init_c_vel**2))

        self.r1_0 = np.array([self.r, 0.0, 0.0])
        self.r2_0 = np.array([0.0, 0.0, 0.0])
        self.v1_0 = np.array([0.0, self.init_e_vel, 0.0])
        self.v2_0 = np.array([0.0, -self.init_e_vel, 0.0])

        self.initial_value = np.hstack((self.r1_0, self.r2_0,
                                        self.v1_0, self.v2_0))
        self.y_noRot = np.ones((tl, self.__eq_amount))

        self.noRot_result = None
        self.rot_result = None
        self.no_rot_energy = None
        self.no_rot_energy_0 = None
        self.energy = None
        self.energy_0 = None

    def __two_body_system(self, eq_index, tmp_array, mu):

        xi = tmp_array[3] - tmp_array[0]
        yi = tmp_array[4] - tmp_array[1]
        zi = tmp_array[5] - tmp_array[2]

        r = np.sqrt(xi**2 + yi**2 + zi**2)

        # particle 1 orbit
        if eq_index == 0:
            return tmp_array[6]
        elif eq_index == 1:
            return tmp_array[7]
        elif eq_index == 2:
            return tmp_array[8]

        # particle 2 orbit
        elif eq_index == 3:
            return tmp_array[9]
        elif eq_index == 4:
            return tmp_array[10]
        elif eq_index == 5:
            return tmp_array[11]

        # particle 1 velocity
        elif eq_index == 6:
            return mu * xi / (r ** 3) + self.potential(self.c, xi, r)
        elif eq_index == 7:
            return mu * yi / (r ** 3) + self.potential(self.c, yi, r)
        elif eq_index == 8:
            return mu * zi / (r ** 3) + self.potential(self.c, zi, r)

        # particle 2 velocity
        elif eq_index == 9:
            return -mu * xi / (r ** 3) - self.potential(self.c, xi, r)
        elif eq_index == 10:
            return -mu * yi / (r ** 3) - self.potential(self.c, yi, r)
        elif eq_index == 11:
            return -mu * zi / (r ** 3) - self.potential(self.c, zi, r)

        else:
            return 0

    def rk4_process(self):
        m1, m2 = self.bh_mass['m1'], self.bh_mass['m2']
        mu = self.__G*(m1 + m2)

        k1, k2 = np.ones(self.__eq_amount), np.ones(self.__eq_amount)
        k3, k4 = np.ones(self.__eq_amount), np.ones(self.__eq_amount)
        v1, v2, v3 = np.ones(self.__eq_amount), np.ones(self.__eq_amount), np.ones(self.__eq_amount)
        time_tmp = 1

        for index in range(self.__eq_amount):
            self.y_noRot[0, index] = self.initial_value[index]

        while time_tmp < self.__tl:
            for index in range(self.__eq_amount):
                k1[index] = self.__dt*self.__two_body_system(index, self.y_noRot[time_tmp - 1], mu)
                v1[index] = self.y_noRot[time_tmp - 1][index] + 0.5*k1[index]
            for index in range(self.__eq_amount):
                k2[index] = self.__dt*self.__two_body_system(index, v1, mu)
                v2[index] = self.y_noRot[time_tmp - 1][index] + 0.5*k2[index]
            for index in range(self.__eq_amount):
                k3[index] = self.__dt*self.__two_body_system(index, v2, mu)
                v3[index] = self.y_noRot[time_tmp - 1][index] + k3[index]
            for index in range(self.__eq_amount):
                k4[index] = self.__dt*self.__two_body_system(index, v3, mu)

            for index in range(self.__eq_amount):
                estimate_term = (k1[index] + 2*k2[index] + 2*k3[index] + k4[index]) / 6.0
                self.y_noRot[time_tmp][index] = self.y_noRot[time_tmp - 1][index] + estimate_term

            time_tmp += 1
        # print("rk4 process done !")

    def no_rot_result(self):
        m1, m2 = self.bh_mass['m1'], self.bh_mass['m2']
        self.rk4_process()
        r1, r2 = self.y_noRot[:, :3], self.y_noRot[:, 3:6]
        v1, v2 = self.y_noRot[:, 6:9], self.y_noRot[:, 9:]

        barycenter = (m1*r1 + m2*r2) / (m1 + m2)

        r1_relation_com = r1 - barycenter
        r2_relation_com = r2 - barycenter
        v1_relation_com = v1
        v2_relation_com = v2
        no_rot_dict = {'p1_orbit': r1_relation_com,
                       'p2_orbit': r2_relation_com,
                       'p1_velocity': v1_relation_com,
                       'p2_velocity': v2_relation_com}

        self.noRot_result = no_rot_dict

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
        no_rot_dict = self.noRot_result

        r1 = no_rot_dict['p1_orbit']
        r2 = no_rot_dict['p2_orbit']
        v1 = no_rot_dict['p1_velocity']
        v2 = no_rot_dict['p2_velocity']

        self.rotation_mat()
        rot_r1 = np.dot(self.rot_mat, r1.T).T
        rot_r2 = np.dot(self.rot_mat, r2.T).T
        rot_v1 = np.dot(self.rot_mat, v1.T).T
        rot_v2 = np.dot(self.rot_mat, v2.T).T

        self.rot_result = {'p1_orbit': rot_r1,
                           'p2_orbit': rot_r2,
                           'p1_velocity': rot_v1,
                           'p2_velocity': rot_v2}

    def cal_total_energy(self):
        no_rot_data = self.noRot_result
        rot_data = self.rot_result
        mass_sum = self.bh_mass['m1'] + self.bh_mass['m2']

        nx1, ny1, nz1 = no_rot_data['p1_orbit'][:, 0], no_rot_data['p1_orbit'][:, 1], no_rot_data['p1_orbit'][:, 2]
        nx2, ny2, nz2 = no_rot_data['p2_orbit'][:, 0], no_rot_data['p2_orbit'][:, 1], no_rot_data['p2_orbit'][:, 2]
        x1, y1, z1 = rot_data['p1_orbit'][:, 0], rot_data['p1_orbit'][:, 1], rot_data['p1_orbit'][:, 2]
        x2, y2, z2 = rot_data['p2_orbit'][:, 0], rot_data['p2_orbit'][:, 1], rot_data['p2_orbit'][:, 2]

        nr1 = np.sqrt(nx1**2 + ny1**2 + nz1**2)
        nr2 = np.sqrt(nx2**2 + ny2**2 + nz2**2)
        r1 = np.sqrt(x1**2 + y1**2 + z1**2)
        r2 = np.sqrt(x2**2 + y2**2 + z2**2)

        n_xv1 = no_rot_data['p1_velocity'][0:, 0]
        n_yv1 = no_rot_data['p1_velocity'][0:, 1]
        n_zv1 = no_rot_data['p1_velocity'][0:, 2]
        n_xv2 = no_rot_data['p2_velocity'][0:, 0]
        n_yv2 = no_rot_data['p2_velocity'][0:, 1]
        n_zv2 = no_rot_data['p2_velocity'][0:, 2]

        xv1 = rot_data['p1_velocity'][0:, 0]
        yv1 = rot_data['p1_velocity'][0:, 1]
        zv1 = rot_data['p1_velocity'][0:, 2]
        xv2 = rot_data['p2_velocity'][0:, 0]
        yv2 = rot_data['p2_velocity'][0:, 1]
        zv2 = rot_data['p2_velocity'][0:, 2]

        no_rotE1 = 0.5*(n_xv1**2 + n_yv1**2 + n_zv1**2) - (mass_sum/nr1) - self.c*(np.arctan(nr1)/nr1)
        no_rotE2 = 0.5*(n_xv2**2 + n_yv2**2 + n_zv2**2) - (mass_sum/nr2) - self.c*(np.arctan(nr2)/nr2)
        no_rot_total_E = no_rotE1 + no_rotE2

        rotE1 = 0.5*(xv1**2 + yv1**2 + zv1**2) - (mass_sum/r1) - self.c*(np.arctan(r1)/r1)
        rotE2 = 0.5*(xv2**2 + yv2**2 + zv2**2) - (mass_sum/r2) - self.c*(np.arctan(r2)/r2)
        rot_total_E = rotE1 + rotE2

        self.no_rot_energy = no_rot_total_E
        self.no_rot_energy_0 = no_rot_total_E[0]
        self.energy = rot_total_E
        self.energy_0 = rot_total_E[0]

    def run(self):
        self.no_rot_result()
        self.rotation_data()
        self.cal_total_energy()

        all_result = {'no_rot_data': self.noRot_result,
                      'rot_data': self.rot_result,
                      'no_rot_energy': self.no_rot_energy,
                      'no_rot_energy0': self.no_rot_energy_0,
                      'rot_energy': self.energy,
                      'rot_energy0': self.energy_0}

        return all_result

    def test(self):
        self.no_rot_result()
        self.rotation_data()
        self.cal_total_energy()

        all_result = {'no_rot_data': self.noRot_result,
                      'rot_data': self.rot_result,
                      'no_rot_energy': self.no_rot_energy,
                      'no_rot_energy0': self.no_rot_energy_0,
                      'rot_energy': self.energy,
                      'rot_energy0': self.energy_0}

        return all_result
