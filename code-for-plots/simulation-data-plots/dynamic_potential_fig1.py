# Plot the shape of the gravitational potential between Andromeda and the Milky Way
#
# First import some modules

from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
from matplotlib import ticker, cm, colors
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from scipy.interpolate import interp1d


# Formula for total potential:
#
#
# V(r,t) = 4 * np.pi * G * (rho_0 * R_s1 ** 3 * (mag_r / (R_s1 + mag_r) - np.log((R_s1 + mag_r) / R_s1)) / mag_r ** 2)
# + rho_0 * R_s2 ** 3 * (mag_andr / (R_s2 + mag_andr) - np.log((R_s2 + mag_andr) / R_s2)) mag_andr ** 2) -
# M1 * mag_r ** 2 / (4 * np.pi * (mag_r ** 2 + a ** 2) ** 1.5) - M2 * mag_andr ** 2 / (4 * np.pi * (mag_andr ** 2 + a ** 2) ** 1.5)))


# Now set the constants, variables and initial conditions

#
'''CONSTANT PARAMETERS:'''
#

# for R_max = 200kpc
rho_0 = 1.948719425e7  # local dark matter density [M_o/kpc^3]

G = 4.49948902 * 10 ** (-12)  # rescaled grav. constant [kpc^3/(Myr^2 * M_o)]

c = 10.  # scale parameter []

R_s1 = 20.  # Milky Way scale radius [kpc]

R_s2 = 20.  # Andromeda scale radius [kpc]

# Plummer radius of galaxies for Plummer term in differential eq., considers only baryonic matter [kpc]
a = 4

# total masses of each galaxy [M_o]
M_MW = 0.8 * 10 ** 12  # mass of the Milky Way
M_AND = 0.8 * 10 ** 12  # mass of Andromeda
# M_MW = 0.615 * 10 ** 12  # mass of the Milky Way
# M_AND = 1.1388 * 10 ** 12  # mass of Andromeda

R_max = 1039.08142  # maximum (radial) distance between MW and AND [kpc]

# dark matter density off halo at scale radius [M_o/kpc^3] (for R_s = 20kpc and M = 1/2 * 5/6 * M_tot)
# rho_0 = 1.1722009 * 10 ** 7

# M_tot = 4.21066612 * 10 ** 12  # total mass of both galaxies [M_o]

# baryonic masses of the two galaxies
M1 = M_MW / 6
M2 = M_AND / 6

# rho_0 if the MW and Andromeda masses are different
rho_0_MW = (5 / 6) * (M_MW / (4 * np.pi * (R_s1 ** 3))) * \
    (np.log(1 + c) - (c / (1 + c))) ** (-1)

rho_0_AND = (5 / 6) * (M_AND / (4 * np.pi * (R_s2 ** 3))) * \
    (np.log(1 + c) - (c / (1 + c))) ** (-1)


# define arrays for the x and y axes
x_axis = np.linspace(-700., R_max + 700., 3000)

# print(x_axis[0:20])

y_axis = np.linspace(-850., 850, 3000)

X, Y = np.meshgrid(x_axis, y_axis)


# print(radiusfunc(t[0]))

# plt.plot(radiusfunc(t))
# plt.show()


# Defining the gravitational potential
def potential(rx, ry, current_separation, G, R_s1, R_s2, M1, M2):
    # magnitude of the distance vector to Milky Way
    mag_r = (rx ** 2 + ry ** 2) ** 0.5

    # magnitude of distance vector to Andromeda
    mag_andr = ((rx - current_separation) ** 2 +
                (ry - 0) ** 2) ** 0.5

    V = (-4) * np.pi * G * (rho_0_MW * R_s1 ** 3 * (np.log(R_s1 + mag_r / R_s1) - mag_r / (R_s1 + mag_r)) / mag_r \
                          + rho_0_AND * R_s2 ** 3 * (np.log(R_s2 + mag_andr / R_s2) - mag_andr / (R_s2 + mag_andr)) / mag_andr \
                          + M1 * mag_r ** 2 / (4 * np.pi * (mag_r ** 2 + a ** 2) ** 1.5) \
                          + M2 * mag_andr ** 2 / (4 * np.pi * (mag_andr ** 2 + a ** 2) ** 1.5))
    
    # V(r,t) = 4 * np.pi * G * (rho_0 * R_s1 ** 3 * (mag_r / (R_s1 + mag_r) - np.log((R_s1 + mag_r) / R_s1)) / mag_r ** 2)
    #                         + rho_0 * R_s2 ** 3 * (mag_andr / (R_s2 + mag_andr) - np.log((R_s2 + mag_andr) / R_s2)) / mag_andr ** 2) -
    # M1 * mag_r ** 2 / (4 * np.pi * (mag_r ** 2 + a ** 2) ** 1.5) - M2 * mag_andr ** 2 / (4 * np.pi * (mag_andr ** 2 + a ** 2) ** 1.5)))
    # print(dy1dt)

    return V


# potential as the z axis value
Z = potential(X, Y, 1000, G, R_s1, R_s2, M1, M2)
print(Z.shape)

print(potential(1000.1, 0, 1000, G, R_s1, R_s2, M1, M2))

# plot contour of the potential
plt.figure(1)
plt.rcParams['font.size'] = 16
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.03)
grav_pot = ax.contour(X, Y, Z, levels = np.linspace(-1., 0, 200), cmap='plasma_r', norm=colors.SymLogNorm(linthresh=0.04, linscale=0.03))
ax.set_xlabel('Distance [kpc]', size=16)
ax.set_ylabel('Distance [kpc]', size=16)
ax.set_aspect('equal')
cb = plt.colorbar(grav_pot, cax=cax, label='Gravitational potential [a.u.]', ticks=[0, -0.2, -0.4, -0.6, -0.8, -1])
ax.set_xlim(-400, 1400)
plt.tight_layout()
plt.savefig('potential_contour_same.pdf')


fig = plt.figure(2)
ax = plt.axes()

ax.plot(x_axis, potential(x_axis, 0, 1000, G, R_s1, R_s2, M1, M2))
plt.show()

'''
plt.figure(2)
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, Z, cmap='Reds')
ax.set_xlabel('Distance [kpc]')
ax.set_ylabel('Distance [kpc]')
ax.set_zlabel('Potential [a.u.]')
ax.set_title('Shape of the potential')
'''
