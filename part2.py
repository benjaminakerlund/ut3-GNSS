import numpy as np
import numpy.linalg as npl
import scipy.linalg as scl
import pandas as pd
import matplotlib.pyplot as plt


## -----------------------------------------------------------------------------
# Data
work_dir = "data/"

c = 299792458.0 # speed of light in vaccum [m/s]
re = 6378136.0 # Earth radius [m]

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
b_sat = pd.read_csv(work_dir+'bsat.csv',delimiter=';')

sat_names = x_sat.columns.to_list()[2:-1]
epochs = x_sat.values[:,1]
x_sat = x_sat.values[:,2:-1].T
y_sat = y_sat.values[:,2:-1].T
z_sat = z_sat.values[:,2:-1].T
b_sat = b_sat.values[:,2:-1].T

n_epochs = epochs.shape[0]
n_sat = x_sat.shape[0]

measures = pd.read_csv(work_dir+'Csimu.csv',delimiter=';').values[:,2:-1].T

## -----------------------------------------------------------------------------
# Position and clock solution computation
'''Copied from part 1...'''
def cosine_matrix(x_sta, y_sta, z_sta, x_sat, y_sat, z_sat):
    n_sat = x_sat.shape[0]

    r = np.sqrt((x_sta-x_sat)**2 + (y_sta - y_sat)**2 + (z_sta - z_sat)**2)
    H = np.vstack(((x_sta - x_sat)/r, (y_sta-y_sat)/r, (z_sta-z_sat)/r, np.ones(n_sat))).T

    return H

sol_lse = np.zeros((n_epochs,4))
residuals = np.zeros((n_sat,n_epochs))
dop = np.zeros((n_epochs,4))

# If we want to consider only a subset of the satellites, change the value of idx_sat
idx_sat = np.arange(n_sat)
# idx_sat = np.array([0,3,4,5,6,7,8])
# idx_sat = np.array([1,2,3,4,5,6,7])
# idx_sat = np.array([0,1,2,4,5,6,7,8])

# Initialisation of the iterative least square algorithm at each epoch
p0 = np.array([0,0,0,0])

# Temporal loop
for it in range(n_epochs):

    # Initialisation
    meas_it = measures[idx_sat, it] # measures at epoch it
    pk = p0.copy() # initialisation of the least_squares
    dp_norm = np.inf # the norm of the residuals
    n_iter = 0 # the current iteration
    iter_max = 20 #the maximum number of iterations
    dp_thr = 1e-6 # the threshold for the residuals

    # iteration
    while  n_iter < iter_max and dp_norm > dp_thr:
        # update iteration number
        n_iter += 1

        # compute measurement residuals
        mes_residuals = meas_it - (np.sqrt((x_sat[idx_sat,it]-pk[0])**2 + (y_sat[idx_sat,it]-pk[1])**2 \
                                        + (z_sat[idx_sat,it]-pk[2])**2) + pk[3] - b_sat[idx_sat,it])

        # direction cosine matrix
        H = cosine_matrix(x_ref, y_ref, z_ref, x_sat[idx_sat, it]-pk[0], y_sat[idx_sat, it]-pk[1], z_sat[idx_sat, it]-pk[2])

        # Solve for the increment of the unknown vector
        dp = scl.inv(H.T @ H) @ H.T @ mes_residuals

        # update the unknown vector
        pk = pk + dp

        dp_norm = scl.norm(dp)

    # message if the iteration numbe reached the limit
    if n_iter == iter_max:
        print("Warning: LSE algorithm has not converged")

    # Compute position, direction cosine matrix and measurement residual with the
    # final position
    H_final = cosine_matrix(x_ref, y_ref, z_ref, x_sat[idx_sat, it]-pk[0], y_sat[idx_sat, it]-pk[1], z_sat[idx_sat, it]-pk[2])
    mes_residuals_final = meas_it - (np.sqrt((x_sat[idx_sat, it] - pk[0]) ** 2 + (y_sat[idx_sat, it] - pk[1]) ** 2 \
                                       + (z_sat[idx_sat, it] - pk[2]) ** 2) + pk[3] - b_sat[idx_sat, it])

    # Compute the DOPs of the solution
    # R_enu =
    # pdop =
    # tdop =
    # vdop =
    # hdop =

    # Store the solution
    # dop[it] = np.array([pdop, tdop, hdop, vdop])
    sol_lse[it] = pk.copy()
    residuals[idx_sat, it] = mes_residuals_final

# Compute the solution (position/clock) error
err_x = sol_lse[:,0] - x_ref
err_y = sol_lse[:,1] - y_ref
err_z = sol_lse[:,2] - z_ref
err_3d = np.sqrt(err_x**2 + err_y**2 + err_z**2)
err_b = sol_lse[:,3] - b0 - b1*(epochs-epochs[0])

##
# plot the solution

plt.figure()
for i_sat in range(n_sat):
    plt.plot(epochs-epochs[0],residuals[i_sat], '.', markersize=2, label=sat_names[i_sat])
plt.xlabel("Time")
plt.ylabel("Measurement residuals [m]")
plt.grid(True)
plt.legend()
plt.draw()

# plot position error
plt.figure()
plt.plot(epochs-epochs[0],err_x,'.', markersize=2, label="error x")
plt.plot(epochs-epochs[0],err_y,'.', markersize=2, label="error y")
plt.plot(epochs-epochs[0],err_z,'.', markersize=2, label="error z")
plt.plot(epochs-epochs[0],err_3d,'.', markersize=2, label="error 3d")
plt.xlabel("Time")
plt.ylabel("Position estimation error [m]")
plt.grid(True)
plt.legend()
plt.draw()

# plot the clock offset error
plt.figure()
plt.plot(epochs-epochs[0],err_b,'.', markersize=2)
plt.xlabel("Time")
plt.ylabel("Clock offset error [m]")
plt.grid(True)
plt.draw()

plt.figure()
plt.plot(epochs-epochs[0],sol_lse[:,3]/c*1e6,'.', markersize=2)
plt.xlabel("Time")
plt.ylabel("Estimated Clock [m]")
plt.grid(True)
plt.draw()


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


