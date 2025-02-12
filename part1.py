import math

import numpy as np
import numpy.linalg as npl
import pandas as pd
import matplotlib.pyplot as plt
import scipy.linalg as scl

# imports to make run on Bena's computer...
import matplotlib
matplotlib.use('QtAgg')


## -----------------------------------------------------------------------------
# Data
work_dir = "data/"

c = 299792458.0     # speed of light in vaccum [m/s]
re = 6378136.0      # Earth radius [m]

# reference position (from TLSG00FRA_R_20203400000_15M_01S_MO.rnx)
x_ref = 4628684.4674
y_ref = 119997.1175
z_ref = 4372110.0327
pos_ref = np.array([x_ref, y_ref, z_ref])

# receiver clock b_ref = b0 + b1*t
b0 = c*50e-6  # [m]   Receiver clock bias
b1 = c*0.2e-9  # [m/s] Receiver clock bias drift

# read the data files
x_sat = pd.read_csv(work_dir+'xsat.csv',delimiter=';')
y_sat = pd.read_csv(work_dir+'ysat.csv',delimiter=';')
z_sat = pd.read_csv(work_dir+'zsat.csv',delimiter=';')

sat_names = x_sat.columns.to_list()[2:-1]
epochs = x_sat.values[:,1]
x_sat = x_sat.values[:,2:-1].T
y_sat = y_sat.values[:,2:-1].T
z_sat = z_sat.values[:,2:-1].T

n_epochs = epochs.shape[0]
n_sat = x_sat.shape[0]

measures = pd.read_csv(work_dir+'Csimu.csv',delimiter=';').values[:,2:-1].T

## -----------------------------------------------------------------------------
# Compute the rotation matrix from the ECEF reference frame to the ENU reference
# frame of a given station (with cartesian coordinates in ECEF frame).

def to_sta_enu_rotation_matrix(x_sta, y_sta, z_sta):
    """
    Computation of the rotation matrix from the ECEF reference frame to the
    station local ENU reference frame

    @param x_sta: x station (ECEF)
    @param y_sta: y station (ECEF)
    @param z_sta: z station (ECEF)
    @return: R the rotation matrix
    """

    # Calculate lat and long for the station
        # station coordinates ECEF
    l_sta = np.arctan( y_sta / x_sta )
    phi_sta = np.arctan2( z_sta , np.sqrt(x_sta**2 + y_sta**2) )
    # atan2 or arctan2 needs a "," delimiter and does the division automatically...


    # Computing the matrix
    R = [
            [-np.sin(l_sta), np.cos(l_sta), 0],
            [-np.sin(phi_sta)*np.cos(l_sta), -np.sin(phi_sta)*np.sin(l_sta), np.cos(phi_sta)],
            [np.cos(phi_sta)*np.cos(l_sta), np.cos(phi_sta)*np.sin(l_sta), np.sin(phi_sta)]
    ]

    return R



## -----------------------------------------------------------------------------
# Compute the azimuth and elevation of the satellites

def to_sat_az_el_d(pos_ref, x_sat, y_sat, z_sat):
    """
    Computation of the satellite azimuth, elevation and range from a station
    reference position on Earth.

    The reference station coordinates and the satellites coordinates have to be
    given in the ECEF reference frame.

    @param pos_ref: station reference position numpy.array(shape=(3,))
    @param x_sat: satelite x position (numpy.ndarray(N,))
    @param y_sat: satelite y position (numpy.ndarray(N,))
    @param z_sat: satelite z position (numpy.ndarray(N,))
    @return: az, zl, d the azimuth, elevation, distance from the station to the
    satellite
    """
    R = to_sta_enu_rotation_matrix(pos_ref[0], pos_ref[1], pos_ref[2])

    [delta_E_sta_sat, delta_N_sta_sat, delta_U_sta_sat] = np.matmul(
        R, [x_sat - pos_ref[0], y_sat - pos_ref[1], z_sat - pos_ref[2]]
    )

    d_sta_sat = np.sqrt(delta_E_sta_sat**2 + delta_N_sta_sat**2 + delta_U_sta_sat**2)
    el = np.arcsin( delta_U_sta_sat / d_sta_sat )
    az = np.arctan2( delta_E_sta_sat , delta_N_sta_sat)

    return az, el, d_sta_sat

## -----------------------------------------------------------------------------
# skyplot of the satellites


fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)

for i_sat in range(n_sat):
    '''az = 0  #TODO Change
    el = 0  #TODO Change
    d = 0   #TODO Change
    '''
    [az, el, d] = to_sat_az_el_d(pos_ref, x_sat[i_sat], y_sat[i_sat], z_sat[i_sat])
    ax.plot(az, 90-el*180/np.pi,'.', markersize=3, label=sat_names[i_sat])

ax.grid(True)
ax.set_yticks(range(0, 100, 10))  # Define the yticks
yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
ax.set_yticklabels(yLabel)
ax.set_ylim(0,90)
ax.set_title("Skyplot")
ax.legend()
plt.draw()
plt.show()

## -----------------------------------------------------------------------------
# Direction cosine matrix

def cosine_matrix(x_sta, y_sta, z_sta, x_sat, y_sat, z_sat):
    """
    Cosine direction matrix
    The ith line of H is equal to:
    [(x_sta(i)-x_sat)/r_i, (y_sta(i)-y_sat)/r_i, (z_sta(i)-z_sat)/r_i, 1]

    @param x_sta: double
    @param y_sta: double
    @param z_sta: double
    @param x_sat: ndarray(nsat,)
    @param y_sat: ndarray(nsat,)
    @param z_sat: ndarray(nsat,)
    @return: H, ndarray, the matrix
    """
    n_sat = x_sat.shape[0]

    r = np.sqrt((x_sta-x_sat)**2 + (y_sta - y_sat)**2 + (z_sta - z_sat)**2)
    H = np.vstack(((x_sta - x_sat)/r, (y_sta-y_sat)/r, (z_sta-z_sat)/r, np.ones(n_sat))).T

    return H


## -----------------------------------------------------------------------------
# DOP computation at reference position


R_ecef2enu = to_sta_enu_rotation_matrix(x_ref, y_ref, z_ref)

dop = np.zeros((n_epochs,4))
# If we want to consider only a subset of the satellites, change the value of idx_sat
# idx_sat = np.arange(n_sat)
# idx_sat = np.array([0,3,4,5,6,7,8])
idx_sat = np.array([1,2,3,4,5,6,7])
idx_sat = np.array([0,1,2,4,5,6,7,8])

for it in range(n_epochs):

    Rot_mat = scl.block_diag(R_ecef2enu, 1)

    H = cosine_matrix(x_ref, y_ref, z_ref, x_sat[idx_sat, it], y_sat[idx_sat, it], z_sat[idx_sat, it])
    M_xyz = npl.inv(H.T @ H)
    pdop = np.sqrt(np.sum(np.diag(M_xyz)[:3]))
    tdop = np.sqrt(M_xyz[3, 3])
    M_enu = Rot_mat @ M_xyz @ Rot_mat.T
    vdop = np.sqrt(M_enu[2,2])
    hdop = np.sqrt(np.sum(np.diag(M_enu)[:2]))

    dop[it] = np.array([pdop, tdop, hdop, vdop])


# display the dops
plt.figure()
plt.plot(epochs-epochs[0], dop[:,0], label='PDOP')
plt.plot(epochs-epochs[0], dop[:,1], label='TDOP')
plt.plot(epochs-epochs[0], dop[:,2], label='HDOP')
plt.plot(epochs-epochs[0], dop[:,3], label='VDOP')
plt.xlabel("Time")
plt.ylabel("Dilutions of precision")
plt.grid(True)
plt.legend()
plt.draw()
plt.show()



