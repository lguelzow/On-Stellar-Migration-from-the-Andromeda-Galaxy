import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
from matplotlib import transforms

import numpy as np
import pandas as pd

import astropy
import healpy as hp
import sys
import scipy.stats
from scipy.interpolate import interp2d
from astropy.coordinates import SkyCoord
import astropy.units as u
from mpl_toolkits.mplot3d import Axes3D


# IF THE VELOCITY DIRECTIONS FORM A LINE AT (140°, [-90°,90°]), THE HEADER IN THE DATA FILES NEEDS TO COMMENTED IN


# CONSTANTS

# velocity for which to filter
filter_velocity = 500

# define constant for transformation from mas/yr to km/s
mas_trafo = 4.744213026

# transformation between km/s and kpc/Myr
trafo = 1.023 / 1000  # 1 kpc/Myr = (1000 / 1.023) km/s

# define sun velocity vector in km/s
# sun_vel = [0, 0, 0]
sun_vel = [11.1, 12.24 + 240, 7.25]
LSR = [0, 240, 0]

# x coordinate of Galactic centre in Galactic coordiates
x_GC = 8.178

# scaling factor for size of velocity vectors in 3D plot
scale = 3

# limits of 3D graph axes
limits = 50

limits_plot = 25
# 7

# Andromeda position in galactic coordinates
andromeda_pos = [121.172, -21.573]

andromeda_gal = SkyCoord(andromeda_pos[0], andromeda_pos[1], frame="galactic", unit="deg")


#
'''READ-OUT data_same FILES: '''
#

# readout process nr 1
# Read out the data of all minimum distance_MWs from the file
# comment in whichecer you need


data_same = pd.read_csv('/home/lguelzow/Nextcloud/MA/MA_Paper/HPC_Results/Plotting/present-time-same-mass.txt', delimiter=' ', header=0)
data_half = pd.read_csv('/home/lguelzow/Nextcloud/MA/MA_Paper/HPC_Results/Plotting/present-time-half-mass.txt', delimiter=' ', header=0)

# print(data_same.columns)

print("Reading data into np.arrays...")

# print(data_half)
# print(data_same)

# read out data_same colums into np arrays
distance_MW_same = np.array(data_same['Present_dist_to_MW[kpc]'])

# parallax = np.array(data_same.parallax)

x_vel_same = np.array(data_same['vx[kpc/Myr]'])

y_vel_same = np.array(data_same['vy[kpc/Myr]'])

z_vel_same = np.array(data_same['vz[kpc/Myr]'])

# print(x_vel_same, y_vel_same, z_vel_same)


# velocity magnitude in galactocentric rest frame in km/s
vel_mag_gc_same = np.array([((x_vel_same[x] + sun_vel[0] * trafo) ** 2 + (y_vel_same[x] + sun_vel[1] * trafo) ** 2 + (z_vel_same[x] + sun_vel[2] * trafo) ** 2) ** 0.5 / trafo \
                        for x in range(len(distance_MW_same))])


# transform to velocity directios in Galactocentric coordinates and rest frane
# by adding velocity of the Sun
v_gc_l_same = np.array([np.arctan2(y_vel_same[x] + sun_vel[1] * trafo, x_vel_same[x] + sun_vel[0] * trafo) * (180 / np.pi) for x in range(len(distance_MW_same))])
v_gc_b_same = np.array([np.arcsin((z_vel_same[x] + sun_vel[2] * trafo) / ((x_vel_same[x] + sun_vel[0] * trafo) ** 2 +
                   (y_vel_same[x] + sun_vel[1] * trafo) ** 2 + (z_vel_same[x] + sun_vel[2] * trafo) ** 2) ** 0.5) * (180 / np.pi) for x in range(len(distance_MW_same))])

# import coordinates in astropy funktion "Skycoord"
galactocentric_velocity_same = SkyCoord(v_gc_l_same[:], v_gc_b_same[:], frame="galactic", unit="deg")

# print(v_gc_l_same)





# read out data_half colums into np arrays
distance_MW_half = np.array(data_half['Present_dist_to_MW[kpc]'])

# parallax = np.array(data_half.parallax)

x_vel_half = np.array(data_half['vx[kpc/Myr]'])

y_vel_half = np.array(data_half['vy[kpc/Myr]'])

z_vel_half = np.array(data_half['vz[kpc/Myr]'])


# velocity magnitude in galactocentric rest frame in km/s
vel_mag_gc_half = np.array([((x_vel_half[x] + sun_vel[0] * trafo) ** 2 + (y_vel_half[x] + sun_vel[1] * trafo) ** 2 \
                         + (z_vel_half[x] + sun_vel[2] * trafo) ** 2) ** 0.5 / trafo for x in range(len(distance_MW_half))])

# print(max(vel_mag_gc_same))

# transform to velocity directios in Galactocentric coordinates and rest frane
# by adding velocity of the Sun
v_gc_l_half = np.array([np.arctan2(y_vel_half[x] + sun_vel[1] * trafo, x_vel_half[x] + sun_vel[0] * trafo) * (180 / np.pi) for x in range(len(distance_MW_half))])
v_gc_b_half = np.array([np.arcsin((z_vel_half[x] + sun_vel[2] * trafo) / ((x_vel_half[x] + sun_vel[0] * trafo) ** 2 +
                   (y_vel_half[x] + sun_vel[1] * trafo) ** 2 + (z_vel_half[x] + sun_vel[2] * trafo) ** 2) ** 0.5) * (180 / np.pi) for x in range(len(distance_MW_half))])

# import coordinates in astropy funktion "Skycoord"
galactocentric_velocity_half = SkyCoord(v_gc_l_half[:], v_gc_b_half[:], frame="galactic", unit="deg")

print('Same mass: ', len(galactocentric_velocity_same))
print('Half mass: ', len(galactocentric_velocity_half))

# print(galactocentric_velocity_same)
# print(galactocentric_velocity_half)

bins=(30, 15)


# prepare input for contour plot
Z_same, x_edges_same, y_edges_same = np.histogram2d(np.deg2rad(v_gc_l_same), np.deg2rad(v_gc_b_same), bins=bins)

Z_half, x_edges_half, y_edges_half = np.histogram2d(np.deg2rad(v_gc_l_half), np.deg2rad(v_gc_b_half), bins=bins)

x_mesh_same = np.convolve(x_edges_same, np.ones(2), 'valid') / 2

y_mesh_same = np.convolve(y_edges_same, np.ones(2), 'valid') / 2

x_mesh_half = np.convolve(x_edges_half, np.ones(2), 'valid') / 2

y_mesh_half = np.convolve(y_edges_half, np.ones(2), 'valid') / 2

# print(Z_same)
# print(x_mesh_same.shape)
# print(y_mesh_same.shape)

interpolation_same = interp2d(x_mesh_same, y_mesh_same, Z_same.T, kind="linear", fill_value=np.nan)

interpolation_half = interp2d(x_mesh_half, y_mesh_half, Z_half.T, kind="linear", fill_value=np.nan)

longitude = np.deg2rad(np.linspace(-180, 180, 360))
latitude = np.deg2rad(np.linspace(-90, 90, 180))

interp_grid_same = interpolation_same(longitude, latitude)

interp_grid_half = interpolation_half(longitude, latitude)

sky_grid = np.meshgrid(longitude, latitude)


#
# Plotting
#

rc('font', **{'size':14})

fig = plt.figure(1)
# make mosaic subplots
ax_dict = fig.subplot_mosaic([["skymap", "c_same"], ["skymap", "c_half"]], \
                             per_subplot_kw={"skymap": {'projection': 'mollweide'}}, width_ratios=[1, 0.05], height_ratios=[1,1])

# make main plot skymap
# ax_dict["skymap"].subplot(projection='mollweide')
# plt.rcParams['font.size'] = 10
# add space for colourbars
divider_same = make_axes_locatable(ax_dict["c_same"])
cax_same = divider_same.append_axes("right", size="300%", pad=0.03)

divider_half = make_axes_locatable(ax_dict["c_half"])
cax_half = divider_half.append_axes("right", size="300%", pad=0.03)
ax_dict["skymap"].grid(True)

# now hide empy axes
ax_dict["c_same"].set_visible(False)
ax_dict["c_half"].set_visible(False)

vel_dir_same = ax_dict["skymap"].contour(*sky_grid, 
                                         interp_grid_same, 
                                         cmap='Blues_r', 
                                         alpha=1, 
                                         linewidths=2,
                                         levels=np.round(np.logspace(0, np.log10(np.nanmax(Z_same) - 32), 3)), 
                                         vmax=np.nanmax(Z_same))

vel_dir_half = ax_dict["skymap"].contour(*sky_grid, 
                                         interp_grid_half, 
                                         cmap='Greys_r', 
                                         alpha=1, 
                                         linewidths=2, 
                                         levels=np.round(np.logspace(0, np.log10(np.nanmax(Z_half) - 55), 3)), 
                                         vmax=np.nanmax(Z_half))

ax_dict["skymap"].scatter(andromeda_gal.l.wrap_at('180d').radian, andromeda_gal.b.radian, s=100, \
            label='Andromeda position', marker='*', color='black')
ax_dict["skymap"].legend(loc='upper right')

# colorbars
cb_same = plt.colorbar(vel_dir_same, cax=cax_same, label='Equal-mass')
cb_half = plt.colorbar(vel_dir_half, cax=cax_half, label='Half-mass')

cb_same.lines[0].set_linewidth(8.0)
cb_half.lines[0].set_linewidth(8.0)

fig.suptitle('Concentration of velocity directions of simulation HVSs at $t_0$ \n Distance to Sun $r_{Sun}<50$ kpc, Gal. coord.') # , fontsize=16)
plt.tight_layout
plt.savefig("Present-vel_directions.pdf")
plt.show()

'''
plt.figure(2)
plt.subplot(111, projection='aitoff')
plt.grid(True)
# generate colourmap
cm = plt.cm.get_cmap('plasma')
# colormap for velocity
vel_cm = vel_mag_gc_same[mask_same]
# plot the points as a scatter plot with velocity cm
star_pos_same = plt.scatter(galactocentric_velocity_same[mask_same].l.wrap_at(
    '180d').radian, galactocentric_velocity_same[mask_same].b.radian, s=15, c=vel_cm, cmap=cm, vmin=350, vmax=700, alpha=0.7)
plt.colorbar(
    star_pos_same, label=r'Velocity magnitude $v_{GC}$ [km/s]', fraction=0.03, pad=0.04)
plt.title('Velocity directions of simulation HVSs at $t_0$ \n Distance to Sun $r_{Sun}<50$ kpc \n Same-mass scenario, Gal. coord.',
          loc='right', fontsize=10)
# plt.legend(loc='upper right')
# plt.savefig("Present-vel_directions.pdf")
plt.show()
'''
