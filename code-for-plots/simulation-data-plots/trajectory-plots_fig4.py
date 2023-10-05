# import modules/libraries
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import sys
import astropy
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.cm as cm
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable


#
'''CONSTANTS AND PARAMETERS: '''
#


# age of the universe [Myrs]
t_0 = 13800

# initial time in simulation [Myrs]
t_10000 = 10000

# length of time interval
time_steps = 10001

# concentration parameter
c = 10

# rescaled grav. constant [kpc^3/(Myr^2 * M_o)]
G = 4.49948902 * 10 ** (-12)

# Milky Way scale radius for NFW dark matter profile [kpc]
R_s1 = 20.

# Andromeda scale radius for NFW dark matter profile [kpc]
R_s2 = 20.  # Andromeda scale radius for NFW dark matter profile [kpc]

# Plummer radius of galaxies for Plummer term in differential eq., considers only baryonic matter [kpc]
a = 4

# total masses of each galaxy [M_o], half mass scenario
# M_MW = 0.615 * 10 ** 12  # mass of the Milky Way
# M_AND = 1.1388 * 10 ** 12  # mass of Andromeda

# total masses of each galaxy [M_o], same mass scenario
M_MW = 0.8 * 10 ** 12  # mass of the Milky Way
M_AND = 0.8 * 10 ** 12  # mass of Andromeda

# baryonic masses of the two galaxies
M1 = M_MW / 6
M2 = M_AND / 6

# rho_0 if the MW and Andromeda masses are different
rho_0_MW = (5 / 6) * (M_MW / (4 * np.pi * (R_s1 ** 3))) * \
    (np.log(1 + c) - (c / (1 + c))) ** (-1)

rho_0_AND = (5 / 6) * (M_AND / (4 * np.pi * (R_s2 ** 3))) * \
    (np.log(1 + c) - (c / (1 + c))) ** (-1)

# position of the Sun in Kartesian coordinate system centered on MW centre
r_Sun = [4.003551, 7.13161944, 0.]

# velocity vector of the sun in the Milky Way rest frame (in Galactic Kartesian coordinates) [km/s]
sun_vel = [11.1, 12.24 + 240, 7.25]

# transformation constant (kpc/Myr -> km/s)
trafo = 1000 / 1.023

# both of these angles at present time:
# angle between (simulation) x-axis and MW galactic coordinate (0,0)-axis [rad]
epsilon = 2.08233

# angle between simulation x-y-plane and MW galactic plane [rad]
phi = -0.4328320437

# rotation matrices for coordinate transformation at the end
# rotation matrix for rotation around z-axis
Rot_z_T = np.array([[np.cos(epsilon), np.sin(epsilon), 0],
                    [np.sin(epsilon) * (-1), np.cos(epsilon), 0], [0, 0, 1]])

# rotation matrix for rotation around x-axis
Rot_x_T = np.array([[1, 0, 0], [0, np.cos(phi), np.sin(phi)],
                    [0, np.sin(phi) * (-1), np.cos(phi)]])


#######################################################################


#
'''VARIABLE PARAMETERS: '''
#

# Read the data of the movement of Andromeda from data file
readout1 = np.loadtxt('M31-trajectory-same-mass.txt', dtype=np.str)
# option for HPC array and the specific filesystem
# readout1 = np.loadtxt('/upb/departments/pc2/users/h/hypvel02/Python_files/new_test/M31-position.txt', dtype=np.str)


# Convert the whole table into an array of strings
# call: readout variable
M31_array = np.array(readout1)

# print(len(M31_array))

# Define the array of time steps in the Andromeda table
t_M31 = [None] * (len(M31_array) - 1)

# Get the time values from the readout variable and convert the strings into floats
for x in range(len(M31_array) - 1):
    t_M31[x] = float(M31_array[x + 1, 0])


# Define arrays for components of Andromeda position vector
# x, y and z components
x_Andr = [None] * (len(M31_array) - 1)

y_Andr = [None] * (len(M31_array) - 1)

z_Andr = [None] * (len(M31_array) - 1)

# x, y and z components of the velocity
vx_Andr = [None] * (len(M31_array) - 1)

vy_Andr = [None] * (len(M31_array) - 1)

vz_Andr = [None] * (len(M31_array) - 1)


# readout Andromeda position from readout variable
for x in range(len(M31_array) - 1):
    # Get the x values from the big table and convert the strings into floats
    x_Andr[x] = float(M31_array[x + 1, 1])

    # Get the y values from the big table and convert the strings into floats
    y_Andr[x] = float(M31_array[x + 1, 2])

    # Get the z values from the big table and convert the strings into floats
    z_Andr[x] = float(M31_array[x + 1, 3])

   # Now do the same for the components of the velocity vector
   # Get the x-direction velocity values from the big table and convert the strings into floats
    vx_Andr[x] = float(M31_array[x + 1, 4])

    # Get the y-direction velocity values from the big table and convert the strings into floats
    vy_Andr[x] = float(M31_array[x + 1, 5])

    # Get the z-direction velocity values from the big table and convert the strings into floats
    vz_Andr[x] = float(M31_array[x + 1, 6])

# print(vx_Andr[0])


# interpolate the position and velocity arrays into continuous functions of time with interp1D
A_xfunc = interp1d(t_M31, x_Andr, bounds_error=False, fill_value="extrapolate")

A_yfunc = interp1d(t_M31, y_Andr, bounds_error=False, fill_value="extrapolate")

A_zfunc = interp1d(t_M31, z_Andr, bounds_error=False, fill_value="extrapolate")


# repeat for the velocities
A_vxfunc = interp1d(t_M31, vx_Andr, bounds_error=False,
                    fill_value="extrapolate")

A_vyfunc = interp1d(t_M31, vy_Andr, bounds_error=False,
                    fill_value="extrapolate")

A_vzfunc = interp1d(t_M31, vz_Andr, bounds_error=False,
                    fill_value="extrapolate")


# plt.plot(A_xfunc(t))
# plt.show()


##################################################################


'''INITIAL CONDITIONS:'''
#

# print(t[0])


# RNG = 0: use freely set initial conditions
# RNG = 1: use randomly generated initial conditions from a separately generated file
RNG = 1

if RNG == 1:
    # Read out the randomly generated initial conditions from file
    readout2 = np.loadtxt('full-trajectories.txt', dtype=np.str)
    # readout2 = np.loadtxt('full-trajectories-weird.txt', dtype=np.str)

    # convert into an array of strings
    ini = np.array(readout2)

    t_ini = [float(ini[x, 3]) for x in range(len(ini))]

    # print(t_ini[0])
    # print(t_ini[1])

    # establish time interval with initial time t_ini = 10000 Myrs
    t = [np.linspace(t_ini[i], t_0, time_steps) for i in range(len(ini))]

    # print(t[0])

    # write the distance values into an array and convert them back to floats
    r0 = [[float(ini[x, 12]) + A_xfunc(t_ini[x]), float(ini[x, 13]) + A_yfunc(t_ini[x]),
           float(ini[x, 14]) + A_zfunc(t_ini[x])] for x in range(len(ini))]

    # print(r0[0])
    # print(r0[1])

    # write the distance values into an array and convert them back to floats
    # r0 = [[float(ini[x, 12]), float(ini[x, 13]),
    #        float(ini[x, 14])] for x in range(len(ini))]

    # print(r0[31])

    # write the velocity values into an array and convert them back to floats
    v0 = [[float(ini[x, 16]) + A_vxfunc(t_ini[x]), float(ini[x, 17]) + A_vyfunc(
        t_ini[x]), float(ini[x, 18]) + A_vzfunc(t_ini[x])] for x in range(len(ini))]

    # write the velocity values into an array and convert them back to floats
    # v0 = [[float(ini[x, 16]), float(ini[x, 17]), float(ini[x, 18])]
    #       for x in range(len(ini))]

    # print(v0[0])
    # print(v0[1])


# scenario for time reversed tracking of HVS in the Milky Way
elif RNG == 0:

    t = np.linspace(t_10000, t_0, time_steps)

    # reversing the time array
    # t_reverse = t[::-1]

    # condition for time to run in reverse
    # comment out if time is supposed to run forwards
    # t = t_reverse

    # use initial condiitions that one can set freely

    # test particle initial position vector to MW centre origin [kpc]
    # r0 = [-5, 0, 0]
    # EG No. 1:
    r0 = [-0.8777564335619109, 5.84152931088274, 7.684609420515938]
    # EG No. 2: r0 = [-2.7460940898620056, 2.3630026025870493, 7.989815745991607]
    # EG No. 3: r0 = [-2.813644081223114, 1.8765290327773057, 4.522083627634269]
    # EG No. 4: r0 = [-9.722413160479924, -2.9936399283928185, 4.591453310239966]
    # EG No. 5: r0 = [-1.0513529910905381, 13.849035588456247, 1.9055613431065033]
    # EG No. 6: r0 = [-3.814947619757579, -0.222517170279453, -5.725299185038648]

    # print(r0)

    # test particle initial velocity vector [kpc/Myr] (can set this freely: ~1000 km/s = (-1.022689977))
    # v0 = [-1.3, 0, 0]
    # EG No. 1:
    v0 = [316.8841280491082 / trafo, -633.5956500991813 /
          trafo, -136.96963406455995 / trafo]
    # EG No. 2: v0 = [588.8448862065513 / trafo, 21.2343831426193 / trafo, -354.3736643688576 / trafo]
    # EG No. 3: v0 = [419.77778555775143 / trafo, -278.75384845443006 / trafo, -447.9201050933998 / trafo]
    # EG No. 4: v0 = [485.37016345394875 / trafo, -261.53288122649144 / trafo, -282.3841232991549 / trafo]
    # EG No. 5: v0 = [635.6531002257959 / trafo, -80.51724561470961 / trafo, -2.3797157229426986 / trafo]
    # EG No. 6: v0 = [431.51292895904146 / trafo, -319.1602964242984 / trafo, 346.13085313598384 / trafo]

    # transformation for ICs to be in the simulation coordinate frame
    r0 = np.dot(Rot_x_T, r0)
    r0 = np.dot(Rot_z_T, r0)

    v0 = np.dot(Rot_x_T, v0)
    v0 = np.dot(Rot_z_T, v0)

    # print(v0)


# initial vector [(kpc, kpc/Myr)]
y0 = [[r0[i][0], v0[i][0], r0[i][1], v0[i][1], r0[i][2], v0[i][2]]
      for i in range(len(ini))]

# print(y0[0])


def diff_eq(y0, t, G, R_s1, R_s2, M1, M2):
    # Unpack the vector y = [rx, vx, ry, vy, rz, vz]
    rx, vx, ry, vy, rz, vz = y0

    # magnitude of distance and velocity vector to Milky Way
    mag_r = (rx ** 2 + ry ** 2 + rz ** 2) ** 0.5
    # mag_v = (vx ** 2 + vy ** 2 + vz ** 2) ** 0.5

    # magnitude of distance vector to Andromeda
    mag_andr = ((rx - A_xfunc(t)) ** 2 +
                (ry - A_yfunc(t)) ** 2 + (rz - A_zfunc(t)) ** 2) ** 0.5

    # calculate differential equations [v_i, a_i]
    dy1dt = [vx, 4 * np.pi * G * (rho_0_MW * R_s1 ** 3 * (mag_r / (R_s1 + mag_r) - np.log((R_s1 + mag_r) / R_s1)) * rx / (mag_r ** 3) + rho_0_AND * R_s2 ** 3 * (mag_andr / (R_s2 + mag_andr) - np.log(
        (R_s2 + mag_andr) / R_s2)) * (rx - A_xfunc(t)) / (mag_andr) ** 3 - M1 * rx / (4 * np.pi * (mag_r ** 2 + a ** 2) ** 1.5) - M2 * (rx - A_xfunc(t)) / (4 * np.pi * (mag_andr ** 2 + a ** 2) ** 1.5))]
    # print(dy1dt)

    dy2dt = [vy, 4 * np.pi * G * (rho_0_MW * R_s1 ** 3 * (mag_r / (R_s1 + mag_r) - np.log((R_s1 + mag_r) / R_s1)) * ry / (mag_r ** 3) + rho_0_AND * R_s2 ** 3 * (mag_andr / (R_s2 + mag_andr) - np.log(
        (R_s2 + mag_andr) / R_s2)) * (ry - A_yfunc(t)) / (mag_andr) ** 3 - M1 * ry / (4 * np.pi * (mag_r ** 2 + a ** 2) ** 1.5) - M2 * (ry - A_yfunc(t)) / (4 * np.pi * (mag_andr ** 2 + a ** 2) ** 1.5))]
    # print(dy2dt)

    dy3dt = [vz, 4 * np.pi * G * (rho_0_MW * R_s1 ** 3 * (mag_r / (R_s1 + mag_r) - np.log((R_s1 + mag_r) / R_s1)) * rz / (mag_r ** 3) + rho_0_AND * R_s2 ** 3 * (mag_andr / (R_s2 + mag_andr) - np.log(
        (R_s2 + mag_andr) / R_s2)) * (rz - A_zfunc(t)) / (mag_andr) ** 3 - M1 * rz / (4 * np.pi * (mag_r ** 2 + a ** 2) ** 1.5) - M2 * (rz - A_zfunc(t)) / (4 * np.pi * (mag_andr ** 2 + a ** 2) ** 1.5))]
    # print(dy3dt)DE

    return dy1dt + dy2dt + dy3dt


# Finally find the solution using scipy.odeint
sol = [odeint(diff_eq, y0[i], t[i], args=(G, R_s1, R_s2, M1, M2))
       for i in range(len(ini))]
# This gives a (len(t), 6) array with the columns corresponding to [rx, vx, ry, vy, rz, vz]


distance = [np.sqrt(sol[i][time_steps - 1, 0] ** 2 + sol[i]
                    [time_steps - 1, 2] ** 2 + sol[i][time_steps - 1, 4] ** 2) for i in range(len(ini))]

# print(distance[0])
# print(distance[1])

#########################################################################################


#
'''PLOTTING:'''
#

minima = 10000 / 1000  # min(t_ini)
maxima = 13000 / 1000  # max(t_ini)

norm = matplotlib.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.get_cmap('plasma'))

# for v in t_ini:
#     print(mapper.to_rgba(v))

# really simple grayscale answer
algebra_list = [(x-minima)/(maxima-minima) for x in t_ini]
# let's compare the mapper and the algebra
mapper_list = [mapper.to_rgba(x / 1000) for x in t_ini]

# print(norm)
# print(sol[31][:, 2])
# print(sol[31][:, 4])

# print(A_xfunc(10441.10867665637))


# Plot position of Milky Way
dummy_x = np.linspace(0., 0.000000001, 10001)
dummy_y = [0] * len(t[0])
dummy_z = [0] * len(t[0])

# for marker of MW centre and Andromeda positions in projection graph
x_point_sgr = 0
x_point_Andr = 774
y_points = 0
z_points = 0


# write all trajectory data into list for plotting
x_line = [sol[i][:, 0] for i in range(len(ini))]
y_line = [sol[i][:, 2] for i in range(len(ini))]
z_line = [sol[i][:, 4] for i in range(len(ini))]

# ax.set_aspect('equal')
'''
plt.figure(1)
ax = plt.axes(projection="3d")
plt.xlabel('x [kpc]')
plt.ylabel('y [kpc]')
ax.set_zlabel('z [kpc]')

# loop to plot all trajectories
for i in range(len(ini)):
    # testparticle trajectory
    ax.plot3D(x_line[0], y_line[0], z_line[0], 'blue')

# plot Milky Way and Andromeda trajectories
# Andromeda trajectory
ax.plot3D(x_Andr, y_Andr, z_Andr, 'black')
# Milky way trajectory/position
ax.plot3D(dummy_x, dummy_y, dummy_z, 'black')

ax.set_xlim(-40, 1060)
ax.set_ylim(-40, 40)
ax.set_zlim(-40, 40)


# condition for all 3 axes to be the same scale
# set scale = 1 for equal scales
scale = 1

if scale == 1:
    max_range = np.array([x_line[0].max() - x_line[0].min(), y_line[0].max() -
                          y_line[0].min(), z_line[0].max() - z_line[0].min()]).max() / 2.

    mid_x = (x_line[0].max() + x_line[0].min()) * 0.5
    mid_y = (y_line[0].max() + y_line[0].min()) * 0.5
    mid_z = (z_line[0].max() + z_line[0].min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
'''

# 2D projection plot with x und y axis
plt.figure(2)
plt.rcParams['font.size'] = 10
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.03)
# plot Milky Way and Andromeda trajectories
# Andromeda trajectory
ax.plot(x_Andr, y_Andr, 'black')
# Milky way trajectory/position
ax.plot(dummy_x, dummy_y, 'black')

for i in range(len(ini)):
    ax.plot(x_line[i], y_line[i], c=mapper_list[i])

# plot a single trajectory alone to add a label
# plt.plot(x_line[0], y_line[0], label='HVS trajectories', c=mapper_list[0])

# plot the MW marker
ax.scatter(x_point_sgr, y_points, label='Milky Way',
            color='black', marker='*', s=85, zorder=100)

# plot the Andromeda marker
ax.scatter(x_point_Andr, y_points, label='Andromeda',
            color='black', marker='s', s=75, zorder=100)

ax.set_xlabel(r'$x$-axis distance [kpc]')
ax.set_ylabel(r'$y$-axis distance [kpc]')
ax.set_xlim(-40, 1100)
ax.set_ylim(-570, 570)
plt.colorbar(
    mapper,cax=cax, label='Time of trajectory start [Gyrs]')
# plt.legend(loc='upper left')

# condition for all 3 axes to be the same scale
# set scale = 1 for equal scales
scale = 1

if scale == 1:
    xy_max_range = np.array([x_line[0].max() - x_line[0].min(), y_line[0].max() -
                             y_line[0].min(), z_line[0].max() - z_line[0].min()]).max() / 2.

    xy_mid_x = (x_line[0].max() + x_line[0].min()) * 0.5
    xy_mid_y = (y_line[0].max() + y_line[0].min()) * 0.5
    xy_mid_z = (z_line[0].max() + z_line[0].min()) * 0.5
    ax.set_xlim(xy_mid_x - xy_max_range - 50, xy_mid_x + xy_max_range + 50)
    ax.set_ylim((xy_mid_y - xy_max_range) * 0.1, (xy_mid_y + xy_max_range) * 0.1)
    # ax.set_zlim(mid_z - max_range, mid_z + max_range)

plt.tight_layout()
plt.savefig('Trajectory-projection-xy.pdf')


# 2D projection plot with x und y axis
plt.figure(3)
plt.rcParams['font.size'] = 10
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.03)
# plot Milky Way and Andromeda trajectories
# Andromeda trajectory
ax.plot(x_Andr, z_Andr, 'black')
# Milky way trajectory/position
ax.plot(dummy_x, dummy_z, 'black')

# plot all the trajectories
plot_list = [ax.plot(x_line[i], z_line[i], c=mapper_list[i])
             for i in range(len(ini))]

# plot a single trajectory alone to add a label
ax.plot(x_line[0], z_line[0], label='HVS trajectories', c=mapper_list[0])

# plot the MW marker
ax.scatter(x_point_sgr, z_points, label='Milky Way',
            color='black', marker='*', s=85, zorder=100)

# plot the Andromeda marker
ax.scatter(x_point_Andr, z_points, label='Andromeda trajectory',
            color='black', marker='s', s=75, zorder=100)

# print(plot_list)

ax.set_xlabel(r'$x$-axis distance [kpc]')
ax.set_ylabel(r'$z$-axis distance [kpc]')
ax.set_xlim(-40, 1100)
ax.set_ylim(-570, 570)
plt.colorbar(
    mapper, label='Time of trajectory start [Gyrs]', cax=cax)
ax.legend(loc='upper right')

if scale == 1:
    xz_max_range = np.array([x_line[0].max() - x_line[0].min(), y_line[0].max() -
                             y_line[0].min(), z_line[0].max() - z_line[0].min()]).max() / 2.

    xz_mid_x = (x_line[0].max() + x_line[0].min()) * 0.5
    xz_mid_y = (y_line[0].max() + y_line[0].min()) * 0.5
    xz_mid_z = (z_line[0].max() + z_line[0].min()) * 0.5
    ax.set_xlim(xz_mid_x - xz_max_range - 50, xz_mid_x + xz_max_range + 50)
    ax.set_ylim((xz_mid_z - xz_max_range) * 0.1, (xz_mid_z + xz_max_range) * 0.1)
    # ax.set_zlim(mid_z - max_range, mid_z + max_range)
plt.tight_layout()
plt.savefig('Trajectory-projection-xz.pdf')
plt.show()
