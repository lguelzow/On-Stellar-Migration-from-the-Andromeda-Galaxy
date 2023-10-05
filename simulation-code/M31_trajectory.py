# This file solves the DEs for two galaxy system
# and returns a table with the position and velocity of Andromeda
# over the chosen time frame
# First import some modules

from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

# Differential equations
# are in the definition around line 90


# Now set the constants, variables and initial conditions

#
'''CONSTANTS AND PARAMETERS:'''
#

t_0 = 13800  # age of the universe [Myr]

G = 4.49948902 * 10 ** (-12)  # rescaled grav. constant [kpc^3/(Myr^2 * M_o)]

# M_tot = 4.211301617 * 10 ** 12  # total mass of both galaxies [M_o]

# individual masses of each galaxy [M_o]
# need to define reduced mass if masses aren't the same
# M1 = 0.8 * 10 ** 12   # Mass of Milky Way
# M2 = 0.8 * 10 ** 12  # Mass of Andromeda

M1 = 0.615 * 10 ** 12  # mass of the Milky Way
M2 = 1.13888888 * 10 ** 12  # mass of Andromeda


#
'''INITIAL CONDITIONS:'''
#

# choose time frame for result table
# create array for the time integration to solve for [Myr]
# always need t_0 as starting point bc of initial conditions, regardless of forwards or backwards
t = np.linspace(13800., 9990., 10001)

# it's important to note here that the x-axis of the simulation is the connection line between the MW centre and M31
# while the y axis is chosen so that both the Sun and M31 are on the xy-plane
# all vectors are in the simulation coordinate frame

# Andromeda position vector at present time [kpc]
r0 = [774.0233265, 0., 0.]

# Andromeda velocity vector at present time [kpc/Myr] ([-109.24999734, -14.47310068, 8.93741788] in km/s)
v0 = [-0.1117627473, -0.014805982, 0.009142978491]


# combine r_0 an v_0  to form initial vector y_0
# initial vector [(kpc, kpc/Myr)]
y0 = [r0[0], v0[0], r0[1], v0[1], r0[2], v0[2]]


#
'''CALCULATIONS:'''
#

# define arrays to read out the results from later
#
# array of distance to Milky Way in every timestep to find out the minimum distance
MW_distance = []

# array for spacial coordinates at minimum distance
position = []

# define array for velocities in each timestep for corresponding velocity
veloc = []

# define array for time value in each step for corresponding time
t_red = []


# Define the vector y = [r, v] which will be our solution
#
# We then return a vector of the derivatives dydt = [r'(t), v'(t)]
#
# which solves the 6 DEs (3 for position, 3 for velocity coordinates)


# define function that contains the DEs and writes the results in every time step into the defined arrays
def diff_eq(y0, t):
    # Unpack the vector y = [rx, vx, ry, vy, rz, vz]
    rx, vx, ry, vy, rz, vz = y0

    # magnitude of distance vector to Milky Way
    mag_r = (rx ** 2 + ry ** 2 + rz ** 2) ** 0.5
    mag_v = (vx ** 2 + vy ** 2 + vz ** 2) ** 0.5

    # print(ry)

    # write the step results into arrays
    # append MW distance to array
    MW_distance.append(mag_r)

    # append spacial coordinates to array
    position.append([rx, ry, rz])

    # append velocity magnitude and components to array
    veloc.append([mag_v, vx, vy, vz])

    # append times to array
    t_red.append(t)

    # differential equations in in component form [v_i, a_i]
    dy1dt = [vx, G * (-1) * (M1 + M2) * rx / (mag_r ** 3)]
    # print(dy1dt)

    dy2dt = [vy, G * (-1) * (M1 + M2) * ry / (mag_r ** 3)]
    # print(dy2dt)

    dy3dt = [vz, G * (-1) * (M1 + M2) * rz / (mag_r ** 3)]
    # print(dy3dt)

    # return DEs for ODEint to solve
    return dy1dt + dy2dt + dy3dt


# Finally find the solution using scipy.ODEint

# ODEint integrates the DEs returned by diff_eq in every timestep
# of the given t-array and gives the solution as a table
sol = odeint(diff_eq, y0, t)
# This gives a (10001, 6) array with the columns correspoding to r(t) and v(t) alternatingly
#
# Then we can print out interesting results from the calculation and finally plot the trajectories at the end


#
'''RESULTS:'''
#

# print out the relevant values on the command line

# minimum distance to MW with coordinates
print("Maximum distance to Milky Way: ")
print(np.amax(MW_distance))
# print(MW_distance[1])

print("at coordinates:")
print(position[np.argmax(MW_distance)])

# time at minimum distance
print("at time t: ")
print(t_red[np.argmax(MW_distance)])

# velocity at minimum distance
print("with velocity vector [magnitude, components]: ")
print(veloc[np.argmax(MW_distance)])


# initial velocity
print(" ")
print("compared to its velocity vector today [magnitude, components]: ")
print(veloc[0])


# since we started at t_0 and went backwards in time we have to
# reverse all the arrays so time flows in the right direction
t_forward = t[::-1]
x_forward = sol[:, 0][::-1]
y_forward = sol[:, 2][::-1]
z_forward = sol[:, 4][::-1]
vx_forward = sol[:, 1][::-1]
vy_forward = sol[:, 3][::-1]
vz_forward = sol[:, 5][::-1]

# save the time + all 3 arrays into textfile in 4 columns
np.savetxt("M31-trajectory.txt", [])
# np.savetxt("M31-trajectory-same-mass.txt", [])
# np.savetxt("M31-trajectory-half-mass.txt", [])

# write table top line
f = open("M31-trajectory.txt", "w")
# f = open("M31-trajectory-same-mass.txt", "w")
# f = open("M31-trajectory-half-mass.txt", "w")
f.write("t_universe" + " " + "x-coordinate" + " " +
        "y-coordinate" + " " + "z-coordinate" + " " + "x-velocity" + " " +
        "y-velocity" + " " + "z-velocity" + "\n")

# write result table into file to be read out by main simulation file
for i in range(len(t)):
    f.write(str(t_forward[i]) + " " + str(x_forward[i]) + " " + str(y_forward[i]) + " " + str(z_forward[i])
            + " " + str(vx_forward[i]) + " " + str(vy_forward[i]) + " " + str(vz_forward[i]) + "\n")


#
'''PLOTTING:'''
#

# optional plot of Andromeda trajectory in Milky Way rest frame
'''
plt.figure(1)
plt.plot(t, sol[:, 0], 'b', label='rx(t)')
plt.xlabel('t [Myr]')
plt.ylabel('rx [kpc]')
plt.legend()
plt.figure(2)
plt.plot(t, sol[:, 1], 'r', label='vx(t)')
plt.xlabel('t [Myr]')
plt.ylabel('vx [kpc/Myr]')
plt.legend()

plt.figure(3)
plt.plot(t, sol[:, 2], 'b', label='ry(t)')
plt.xlabel('t [Myr]')
plt.ylabel('ry [kpc]')
plt.legend()
plt.figure(4)
plt.plot(t, sol[:, 3], 'r', label='vy(t)')
plt.xlabel('t [Myr]')
plt.ylabel('vy [kpc/Myr]')
plt.legend()

plt.figure(5)
plt.plot(t, sol[:, 4], 'b', label='rz(t)')
plt.xlabel('t [Myr]')
plt.ylabel('rz [kpc]')
plt.legend()
plt.figure(6)
plt.plot(t, sol[:, 5], 'r', label='vz(t)')
plt.xlabel('t [Myr]')
plt.ylabel('vz [kpc/Myr]')
plt.legend()


# short term arrays for 3D plots
dummy_x = np.linspace(0., 0.0000000001, 10001)
dummy_y = [0] * len(t)
dummy_z = [0] * len(t)

plt.figure(7)
ax = plt.axes(projection="3d")
plt.xlabel('x [kpc]')
plt.ylabel('y [kpc]')
# plt.zlabel('z [kpc]')

x_line = sol[:, 0]
y_line = sol[:, 2]
z_line = sol[:, 4]

# ax.set_aspect('equal')

# testparticle trajectory
ax.plot3D(x_line, y_line, z_line, 'blue')
# Milky way trajectory/position
ax.plot3D(dummy_x, dummy_y, dummy_z, 'black')


# condition for all 3 axes to be the same scale
# set scale = 1 for equal scales
scale = 1

if scale == 1:
    max_range = np.array([x_line.max() - x_line.min(), y_line.max() -
                          y_line.min(), z_line.max() - z_line.min()]).max() / 2.

    mid_x = (x_line.max() + x_line.min()) * 0.5
    mid_y = (y_line.max() + y_line.min()) * 0.5
    mid_z = (z_line.max() + z_line.min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)


plt.show()
'''
