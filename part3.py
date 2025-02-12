import numpy as np
import numpy.linalg as npl
import scipy.linalg as scl
import pandas as pd
import matplotlib.pyplot as plt


## -----------------------------------------------------------------------------
# Data
work_dir =

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

sat_names = x_sat.columns.to_list()[2:-1]
epochs = x_sat.values[:,1]
x_sat = x_sat.values[:,2:-1].T
y_sat = y_sat.values[:,2:-1].T
z_sat = z_sat.values[:,2:-1].T

n_epochs = epochs.shape[0]
n_sat = x_sat.shape[0]

measures = pd.read_csv(work_dir+'Csimu.csv',delimiter=';').values[:,2:-1].T

## -----------------------------------------------------------------------------
# simulate a faulty satellite
faulty_sat = 1
faulty_measures = measures.copy()

# fault with an offset
faulty_measures[faulty_sat] = faulty_measures[faulty_sat] + 40*(epochs - epochs[0] > 400)

# fault with a drift
# faulty_measures[faulty_sat] = faulty_measures[faulty_sat] + 0.03*(epochs - epochs[0] > 400)*(epochs-epochs[0]-400)

## -----------------------------------------------------------------------------
# Solution of positionning problem with LSE algorithm (with the faulty satellite)


sol_lse_fault = np.zeros((n_epochs,4))
residuals_fault = np.zeros((n_sat,n_epochs))
dop_fault = np.zeros((n_epochs,4))

# If we want to consider only a subset of the satellites, change the value of idx_sat
idx_sat = np.arange(n_sat)

# Initialisation of the iterative least square algorithm at each epoch
p0 =


# Temporal loop
for it in range(n_epochs):

    # Initialisation


    # iteration
    while  n_iter < iter_max and dp_norm > dp_thr:


    # message if the iteration numbe reached the limit
    if n_iter == iter_max:
        print("Warning: LSE algorithm has not converged")

    # Compute position, direction cosine matrix and measurement residual with the
    # final position
    H_final =
    mes_residuals_final =

    # Compute the DOPs of the solution

    # Store the solution
    dop_fault[it] = np.array([pdop, tdop, hdop, vdop])
    sol_lse_fault[it] = pk.copy()
    residuals_fault[idx_sat,it] = mes_residuals_final


# Compute position error
err_x_fault = sol_lse_fault[:,0] - x_ref
err_y_fault = sol_lse_fault[:,1] - y_ref
err_z_fault = sol_lse_fault[:,2] - z_ref
err_3d_fault = np.sqrt(err_x_fault**2 + err_y_fault**2 + err_z_fault**2)
err_b_fault = sol_lse_fault[:,3] - b0 - b1*(epochs-epochs[0])

## -----------------------------------------------------------------------------
# Plotting the analysis of the estimation with faulty measurements
# plot residuals

plt.figure()
for i_sat in range(n_sat):
    plt.plot(epochs-epochs[0],residuals_fault[i_sat],label=sat_names[i_sat])
plt.xlabel("Time")
plt.ylabel("Measurement residuals [m]")
plt.grid(True)
plt.legend()
plt.title("Faulty satellite: "+sat_names[faulty_sat])
plt.draw()

# plot position error
plt.figure()
plt.plot(epochs-epochs[0],err_x_fault,label="error x")
plt.plot(epochs-epochs[0],err_y_fault,label="error y")
plt.plot(epochs-epochs[0],err_z_fault,label="error z")
plt.plot(epochs-epochs[0],err_3d_fault,label="error 3d")
plt.xlabel("Time")
plt.ylabel("Position estimation error [m]")
plt.grid(True)
plt.legend()
plt.title("Faulty satellite: "+sat_names[faulty_sat])
plt.draw()

# plot the clock offset error
plt.figure()
plt.plot(epochs-epochs[0],err_b_fault)
plt.xlabel("Time")
plt.ylabel("Clock offset error [m]")
plt.grid(True)
plt.title("Faulty satellite: "+sat_names[faulty_sat])
plt.draw()


plt.figure()
plt.plot(epochs-epochs[0],sol_lse_fault[:,3]/c*1e6)
plt.xlabel("Time")
plt.ylabel("Estimated Clock [m]")
plt.grid(True)
plt.draw()


# display the dops
plt.figure()
plt.plot(epochs-epochs[0], dop_fault[:,0], label='PDOP')
plt.plot(epochs-epochs[0], dop_fault[:,1], label='TDOP')
plt.plot(epochs-epochs[0], dop_fault[:,2], label='HDOP')
plt.plot(epochs-epochs[0], dop_fault[:,3], label='VDOP')
plt.xlabel("Time")
plt.ylabel("Dilutions of precision")
plt.grid(True)
plt.legend()
plt.draw()


## -----------------------------------------------------------------------------
# Fault detection algorithm to (try to) detect the faulty satellite


sol_lse_fde = []
dop_fde = []
residuals_fde = []
err_x_fde = []
err_y_fde = []
err_z_fde = []
err_3d_fde = []
err_b_fde = []
distance_full_solution = []
# res_distance_full_solution = []

# loop over all the satellite to exclude the satellite i_sat
# this is the i-th fault mode
for i_sat in range(n_sat):

    print(f"Excluding satellite {sat_names[i_sat]}")

    sol_lse_mode = np.zeros((n_epochs,4))
    residuals_mode = np.zeros((n_sat, n_epochs))
    dop_mode = np.zeros((n_epochs, 4))

    p_it = np.array([0.0, 0.0, 0.0, 0.0])

    idx_sat = np.where(np.arange(n_sat) != i_sat)[0]

    # Solve the LSE algorithm with all the measures except the measures of the
    # satellite i_sat

    for it in range(n_epochs):


# Temporal loop
for it in range(n_epochs):

    # Initialisation


    # iteration
    while  n_iter < iter_max and dp_norm > dp_thr:


        # message if the iteration numbe reached the limit
        if n_iter == iter_max:
            print("Warning: LSE algorithm has not converged")

        # Compute position, direction cosine matrix and measurement residual with the
        # final position
        H_final =
        mes_residuals_mode_final =

        # Compute the DOPs of the solution

        # Store the solution
        dop_mode[it] = np.array([pdop, tdop, hdop, vdop])
        sol_lse_fault[it] = pk.copy()
        residuals_mode[idx_sat,it] = mes_residuals_final


    # Compute position error
    err_x_fde.append(sol_lse_fault[:, 0] - x_ref)
    err_y_fde.append(sol_lse_fault[:, 1] - y_ref)
    err_z_fde.append(sol_lse_fault[:, 2] - z_ref)
    err_3d_fde.append(np.sqrt(
        err_x_fault ** 2 + err_y_fault ** 2 + err_z_fault ** 2))
    err_b_fde.append(sol_lse_fault[:, 3] - b0 - b1 * (epochs - epochs[0]))

    sol_lse_fde.append(sol_lse_mode)
    dop_fde.append(dop_mode)
    residuals_fde.append(residuals_mode)

    distance_full_solution.append(npl.norm((sol_lse_mode-sol_lse_fault)[:,:2], axis=1))

## Plotting the analysis of the estimation with the fault detection and exclusion

# Distance between the i-th fault mode solution and the solution with all the satellites
plt.figure()
for i_sat in range(n_sat):
    plt.plot(epochs-epochs[0],distance_full_solution[i_sat],label=sat_names[i_sat]+'-th fault mode')
plt.xlabel("Time")
plt.ylabel("Distance between each faulty mode and the all-in-view solution [m]")
plt.grid(True)
plt.legend()
plt.title("Faulty satellite: "+sat_names[faulty_sat])
plt.draw()

# Comparison of the residuals
plt.figure()
for it in range(3):
    for jt in range(3):
        plt.subplot(3,3,3*it+jt+1)
        for i_sat in range(n_sat):
            if i_sat != 3*it+jt:
                plt.plot(epochs-epochs[0],residuals_fde[3*it+jt][i_sat],label=sat_names[i_sat])
        plt.legend()
        plt.title(sat_names[3*it+jt]+"-th fault mode")
        plt.grid(True)

plt.suptitle("Comparison of the residuals for each fault mode")
plt.draw()

