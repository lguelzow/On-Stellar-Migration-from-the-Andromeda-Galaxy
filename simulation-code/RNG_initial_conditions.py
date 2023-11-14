# This file generates random initial conditions for the HVS we simulate in the main file

import numpy as np
import random as rng
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
# import scipy.integrate as integrate


#
'''CONSTANTS AND PARAMETERS:'''
#

# amount of initial conditions generated
ini_con = 250000

# rescaled grav. constant [kpc^3/(Myr^2 * M_o)]
G = 4.49948902 * 10 ** (-12)

# Andromeda scale radius for NFW dark matter profile [kpc]
R_s = 20.

# concentration parameter
c = 10

# max distance for escaping Andromeda
R_max = 3 * R_s * c

# Plummer radius (core of the cluster) [kpc]
plum_r = 4

# transformation constant for km/s -> kpc/Myr
trafo = 1.023 / 1000

# baryonic mass of Andromeda [M_o]
mass_and = (1 / 6) * 0.8 * 10 ** 12
mass_MW = mass_and

# mass variant for Andromeda double the size of the Milky Way
# mass_and = (1 / 6) * 1.138888888 * 10 ** 12
# mass_MW = 0.5 * mass_and

# sqrt of mass ratio of Milky Way and Andromeda
mass_ratio_sqrt = mass_and / mass_MW

# mass of Sagittarius A*
mass_sgt = 4.3 * 10 ** 6

# characteristic density of Dark Matter halo according to NFW profile
rho_0 = (5 / 6) * (6 * mass_and / (4 * np.pi * (R_s ** 3))) * \
    (np.log(1 + c) - (c / (1 + c))) ** (-1)

# maximum baryon density of cluster, calculated from Plummer model: rho(r = 0) [M_o / kpc^3]
rho_max = (3 * mass_and) / (4 * np.pi * plum_r ** 3) * (1 ** (-5 / 2))

# fit parameter (velocity scale)
# sigma_v = 90
sigma_v = 59.8842261

# value at lower cut-off of velocity distribution (maximum)
star_number = (1 / sigma_v) # * np.exp(-1 / sigma_v)


#
'''CALCULATIONS:'''
#


# define function for Plummer radius depending on density
def radius(a, rho, M):

    r = a * ((((4 * np.pi * a ** 3 * rho)/(3 * M)) ** (-0.4)) - 1) ** 0.5
    return r


# define function for gravitational potential in Andromeda depending on distance to center
# with NFW and Plummer terms
def grav_potential(r):

    # combined Plummer and NFW potential
    U = (-1) * G * (mass_and / (np.sqrt(r ** 2 + plum_r ** 2)) +
                    4 / r * np.pi * rho_0 * R_s ** 3 * np.log(1 + (r / R_s)))

    # Plummer and NFW potential with added SMBH point mass
    # U = (-1) * G * (mass_sgt / r + (mass_and - mass_sgt) / (np.sqrt(r ** 2 + plum_r ** 2)) + 4 / r *
    #                np.pi * rho_0 * R_s ** 3 * np.log(1 + (r / R_s)))

    return U


# define integrand for integral over force to escape Andromeda
def integrand(r):

    # force from derivation of wikipedia potentials
    F = G * mass_and * r / (r ** 2 + plum_r ** 2) ** 1.5 + 4 * np.pi * G * \
        rho_0 * R_s ** 3 / r ** 2 * (np.log((r + R_s) / R_s) - r / (r + R_s))

    return F


# test for this function
# test = integrate.quad(integrand, 15, 10000)

# v_escape = np.sqrt(2 * abs(test[0])) / trafo

# print(test)
# print(v_escape)


# define function for Andromeda's escape velocity at a given radius
def escape_vel(r):

    # determine v_esc directly through the potential
    # v_esc = np.sqrt(8 * np.pi * G * (rho_0 * R_s ** 3 / r * np.log(1 +
    #                                                                r / R_s) + mass_and / (4 * np.pi * np.sqrt(r ** 2 + plum_r ** 2))))
    v_esc = np.sqrt((-2) * grav_potential(r))
    # v_esc = np.sqrt(2 * abs(grav_potential(r) - grav_potential(R_max)))

    # escape velocity with integral over force required
    # v_esc = np.sqrt(2 * abs(integrate.quad(integrand, r, R_max)[0]))

    # transform to km/s for result
    return v_esc / trafo


# print("\nEscape velocity at 5kpc: ", escape_vel(5), "\n")
# print("\nEscape velocity at 15kpc: ", escape_vel(15), "\n")


# define function for generation of initial velocities based on the local escape velocity
def ini_vel_old(v_esc, rng_star):

    v = v_esc - (sigma_v * (np.log(rng_star) + np.log(sigma_v)))
    return v

# inverse sampling the CDF of the exponential distribution
def ini_vel_exp(v_esc):

    v = v_esc - (sigma_v * np.log(1 - rng.random()))
    return v

# alternative function where probability for high velocities decreases like a power law
# inverse sample the CDF
def ini_vel_power_law(v_esc):

    v = v_esc * (1 / (1 - rng.random())) ** (1 / (3.9))
    return v


# list of the values of the gravitational potential of Andromeda
# at given radii[kpc]
pot = [grav_potential(x) for x in np.arange(0.1, 40, 0.1)]

ve = [escape_vel(x) for x in np.arange(0.1, 40, 0.1)]

# print(np.amax(ve))


# list of randomly generated radii with Plummer model[kpc]
r_coord = [radius(plum_r, rng.random() * rho_max, mass_and)
           for i in range(ini_con)]

# r_coord = [15] * ini_con


# isotropic angular coordinates
# generate list of random phi/declination coordinates
phi_coord = [rng.uniform(0, 2 * np.pi) for i in range(ini_con)]

# random, but weighted theta coordinates / inclination angles
theta_coord = [np.arccos(rng.uniform(-1, 1)) for i in range(ini_con)]

# generate list of random send-off times[Myrs]
# vary initial send-off time for different Andromeda masses
# 10000 Myrs for m=1.13 Mo AND 10500 Myrs for m=0.9 Mo
send_off_times = [rng.uniform(10000, 13000) for i in range(ini_con)]
# print(theta_test)


# list of energy necessary to escape dependent on r_coord
# escape_energy = [integrate.quad(integrand, r_coord[i], R_max)[0]
#                  for i in range(ini_con)]

# print(escape_energy[0])

# generate list of escape velocities corresponding to r_coord list
v_esc_list = [escape_vel(r_coord[i]) for i in range(ini_con)]

# print(v_esc_list[0:10])

# generate list of initial velocity magnitudes
# mag_vel = [trafo * ini_vel_old(v_esc_list[i], star_number * rng.random()) for i in range(ini_con)]

# generate list of initial velocity magnitudes
mag_vel = [trafo * ini_vel_power_law(v_esc_list[i]) for i in range(ini_con)]

mag_vel_test = [mag_vel[i] / trafo 
                    for i in range(ini_con)]  # if mag_vel[i] / trafo > 800]

# print(len(mag_vel_test))

# generate isotropic angular coordinates for velocity
# generate list of random phi/declination coordinates
phi_vel = [rng.uniform(0, 2 * np.pi) for i in range(ini_con)]

# random, but weighted theta coordinates / inclination angles
theta_vel = [np.arccos(rng.uniform(-1, 1)) for i in range(ini_con)]


# transform spherical coordinates to Cartesian
# first for position vector components

r1 = [r_coord[i] * np.sin(theta_coord[i]) * np.cos(phi_coord[i])
      for i in range(ini_con)]

r2 = [r_coord[i] * np.sin(theta_coord[i]) * np.sin(phi_coord[i])
      for i in range(ini_con)]

r3 = [r_coord[i] * np.cos(theta_coord[i]) for i in range(ini_con)]

# then repeat for velocity vector components
v1 = [mag_vel[i] * np.sin(theta_vel[i]) * np.cos(phi_vel[i])
      for i in range(ini_con)]

v2 = [mag_vel[i] * np.sin(theta_vel[i]) * np.sin(phi_vel[i])
      for i in range(ini_con)]

v3 = [mag_vel[i] * np.cos(theta_vel[i]) for i in range(ini_con)]


#
'''RESULTS: '''
#

# save all 6 initial conditions into a line in the textfile for main simulation file to read
np.savetxt("HVS-initial-conditions.txt", [])

# write column headers
f = open("HVS-initial-conditions.txt", "w")
f.write("initial_time" + " " + "x-position" + " " + "y-position" + " " + "z-position" + " " +
        "x-velocity" + " " + "y-velocity" + " " + "z-velocity" + " " + "velocity_magnitude" "\n")

# write in results
for i in range(ini_con):
    # write in the cartesian coordinates for distance and velocity in each line
    f.write(str(send_off_times[i]) + " " + str(r1[i]) + " " + str(r2[i]) + " " + str(r3[i]) + " " +
            str(v1[i]) + " " + str(v2[i]) + " " + str(v3[i]) + " " + str(mag_vel[i]) + "\n")


#
'''PLOTTING:'''
#

# optional plots for inital velocity, time and directions

# plot distribution of radii
plt.figure(1)
# plt.hist(test_plot, 30)
plt.hist(r_coord, bins=75)
plt.xlabel('Initial radius [kpc]')
plt.ylabel('% of test stars')
# plt.yscale('log')
# plt.xlim(0, 10)
plt.savefig('initial_rad.pdf')


# plot distribution of initial velocities
plt.figure(2)
# plt.hist(test_plot, 30)
plt.hist(mag_vel_test, bins=1000)
plt.xlabel('Initial velocity magnitude [km/s]')
plt.ylabel('% of test stars')
# plt.yscale('log')
# plt.xscale('log')
plt.xlim(400, 2500)
plt.savefig('initial_vel_mag.pdf')


# plot of potential level
fig = plt.figure(3)
ax = plt.axes()

ax.plot(np.arange(0.1, 40, 0.1), ve)
plt.xlim(0, 35)
plt.ylim(200, 900)
plt.show()

'''
# kartesian plot of star positions and velocities
X, Y, Z, U, V, W = zip(*test_plot[0:100])

plt.figure(2)
ax = plt.subplot(111, projection='3d')
ax.quiver(X, Y, Z, U, V, W)
ax.set_xlim([-5, 5])
ax.set_ylim([-5, 5])
ax.set_zlim([-5, 5])
ax.set_xlabel('x-axis [kpc]')
ax.set_ylabel('y-axis [kpc]')
ax.set_zlabel('z-axis [kpc]')
plt.show()


plt.hist(t_coord, 80)
plt.xlabel('Distance to center of Andromeda [kpc')
plt.ylabel('# of results')
plt.show()
'''
