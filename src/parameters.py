from turtle import distance
import numpy as np

"""
Mass unit   : 1.2e+12 (Sun mass)
Length unit : kpc
Time unit   : 4.3e+15 (Years)
G = 1

Mass of Galaxy (M87)    : ((4.8~7) * 1e+12).  Choose 5.7e+12 (Sun Mass)
Black Hole (M87)        : ((6.2~7)*1e+9).     Choose 6.4e+6  (Sun Mass)
"""



# radian : 
rot_para = {"omega" : np.pi / 6, 
            "I" : np.pi / 4, 
            "Omega" : np.pi / 6}


init_para = {"black_hole1_mass" : 0.5,  # Sun Mass
             "black_hole2_mass" : 0.5, 
             "gravitational_constant" : 1,
             "distance" : 0.2,   # kpc 
             "eccentricity" : 0.5, 
             "galactic_mass" : (5.7e+12 / 6.4e+9) * 0.5}


mass_exp = {"m1_1" : 0.5, 
            "m1_2" : 0.6, 
            "m1_3" : 0.7, 
            "m1_4" : 0.8, 
            "m1_5" : 0.9, 
            "m2": 0.5}


# distance + epsilon : np.array([r1's epsilon, r2's epsilon])
distance_exp1 = {"r_1" : np.array([0.0, 0.0]), 
                 "r_2" : np.array([0.1, 0.0]), 
                 "r_3" : np.array([0.2, 0.0]), 
                 "r_4" : np.array([0.3, 0.0]), 
                 "r_5" : np.array([0.4, 0.0])}


distance_exp2 = {"r_1" : np.array([0.0, 0.0]), 
                 "r_2" : np.array([0.1, 0.1]), 
                 "r_3" : np.array([0.2, 0.2]), 
                 "r_4" : np.array([0.3, 0.3]), 
                 "r_5" : np.array([0.4, 0.4])}


eccentricity_exp = {"ecc_1" : 0.3, 
                    "ecc_2" : 0.4, 
                    "ecc_3" : 0.5, 
                    "ecc_4" : 0.6, 
                    "ecc_5" : 0.7}



