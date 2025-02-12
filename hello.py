from part1 import *
import numpy as np


# imports to make run on Bena's computer...
import matplotlib
matplotlib.use('QtAgg')


'''Constants
'''
c = 299792458 # [m/s]

'''Receiver position
* receiver is static
* 
'''
x_sta = 4628684.4674    # [m]
y_sta = 119997.1175     # [m]
z_sta = 4372110.0327    # [m]
pos_ref = np.array([x_sta, y_sta, z_sta])



print("Testing rotation matrix function:")
R = to_sta_enu_rotation_matrix(x_sta, y_sta, z_sta)
print(R)

print("Testing azimuth, elevation and distance")
[az, el, d_sta_sat] = to_sat_az_el_d(pos_ref, x_sta, y_sta, z_sta) # using for ref station,
print("Azimuth:" + az + ", Elevation" + el + ", Distance: " + d_sta_sat)