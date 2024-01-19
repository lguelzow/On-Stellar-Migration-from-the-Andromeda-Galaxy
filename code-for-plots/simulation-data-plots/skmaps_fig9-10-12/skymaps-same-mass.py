import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib.ticker import PercentFormatter
from mpl_toolkits import mplot3d
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors
import numpy as np
import astropy
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.optimize import curve_fit


# constant for transformation from kpc/Myr to km/s
trafo = 1000 / 1.023  # 1 kpc/Myr = (1000 / 1.023) km/s

# scaaling factor for 3D plot
scale = 3

# sun velocity vector in galactic rest frame
sun_vel = [11.1, 12.24 + 240, 7.25]

# Andromeda position in galactic coordinates
andromeda_pos = [121.172, -21.573]

andromeda_gal = SkyCoord(
    andromeda_pos[0], andromeda_pos[1], frame="galactic", unit="deg")

#
'''READ-OUT DATA FILES: '''
#

# present distances ----------------------------------------------------------------

# readout process nr 1
# Read out the data of all distances at t_0
readout1 = np.loadtxt('/home/lguelzow/Nextcloud/MA/MA_Paper/Paper_github/On-Stellar-Migration-from-the-Andromeda-Galaxy/simulation-results/power_law_test/combined/present-time-power-law.txt', dtype=np.str)
# readout1 = np.loadtxt('present-time-same-mass.txt', dtype=np.str)

# Convert the whole table into an array of strings
dist1 = np.array(readout1)

# Get the distance values from the read-out table and convert the strings into floats
present_sun_distance = [float(dist1[x, 1])
                        for x in range(len(dist1)) if float(dist1[x, 1]) < 40]

# Get the distance values from the read-out table and convert the strings into floats
present_MW_distance = [float(dist1[x, 0]) for x in range(len(dist1))]

# get present velocity for colormap
present_vel_mag = [trafo * np.sqrt(float(dist1[x, 9]) ** 2 + float(dist1[x, 10])
                                   ** 2 + float(dist1[x, 11]) ** 2) for x in range(len(dist1)) if float(dist1[x, 1]) < 40]

# and again for the MW distance maps
present_vel_mag_MW = [trafo * np.sqrt(float(dist1[x, 9]) ** 2 + float(dist1[x, 10])
                                      ** 2 + float(dist1[x, 11]) ** 2) for x in range(len(dist1))]

# print(np.amin(present_vel_mag))
# print(np.amax(present_vel_mag))

# extract galactic longitude
present_gal_l = [float(dist1[x, 4])
                 for x in range(len(dist1)) if float(dist1[x, 1]) < 40]

# extract galactic latitude
present_gal_b = [float(dist1[x, 5])
                 for x in range(len(dist1)) if float(dist1[x, 1]) < 40]

# now feed the coordinate lists into skycoord arrays
present_gal = SkyCoord(
    present_gal_l[:], present_gal_b[:], frame="galactic", unit="deg")


# do full length lists for trafo to Galactocentric system
# sun distance
present_sun_distance_trafo = [float(dist1[x, 1]) for x in range(len(dist1))]

# extract galactic longitude
present_gal_l_trafo = [float(dist1[x, 4]) for x in range(len(dist1))]

# extract galactic latitude
present_gal_b_trafo = [float(dist1[x, 5]) for x in range(len(dist1))]


# transformation to radial ccordinates around MW centre
present_x = [8.17827 - (present_sun_distance_trafo[i] * np.cos(present_gal_b_trafo[i] * np.pi / 180)
                        * np.cos(present_gal_l_trafo[i] * np.pi / 180)) for i in range(len(dist1))]

present_y = [present_sun_distance_trafo[i] * np.cos(present_gal_b_trafo[i] * np.pi / 180) * np.sin(
    present_gal_l_trafo[i] * np.pi / 180) for i in range(len(dist1))]

present_z = [present_sun_distance_trafo[i] *
             np.sin(present_gal_b_trafo[i] * np.pi / 180) for i in range(len(dist1))]

present_r = [np.sqrt(present_x[i] ** 2 + present_y[i] **
                     2 + present_z[i] ** 2) for i in range(len(dist1))]

# and transform back to spherical for now Galactocentric coordinates
present_galcentre_l = [np.arctan2(present_y[i], present_x[i])
                       * (180 / np.pi) for i in range(len(dist1))]

present_galcentre_b = [np.arcsin(present_z[i] / present_r[i])
                       * (180 / np.pi) for i in range(len(dist1))]

# print(present_r[0:10])
# print(present_MW_distance[0:10])

# now feed the coordinate lists into skycoord arrays
present_galcentre = SkyCoord(
    present_galcentre_l[:], present_galcentre_b[:], frame="galactic", unit="deg")

# print(present_gal[0:10])
# print(present_galcentre_l[0:10])


# minimum distances -----------------------------------------------------------------

# readout process nr 1
# Read out the data of all distances at t_0
# readout1 = np.loadtxt('all-min-distances-to-sun.txt', dtype=np.str)
readout2 = np.loadtxt('/home/lguelzow/Nextcloud/MA/MA_Paper/Paper_github/On-Stellar-Migration-from-the-Andromeda-Galaxy/simulation-results/power_law_test/combined/present-time-power-law.txt', dtype=np.str)

# Convert the whole table into an array of strings
dist2 = np.array(readout2)

print("\n # of present time results:" + str(len(dist1)) + "\n")
print("\n # of minimum distance results:" + str(len(dist2)) + "\n")

# Get the distance values from the read-out table and convert the strings into floats
min_sun_distance = [float(dist2[x, 1])
                    for x in range(len(dist2)) if float(dist2[x, 1]) < 40]

# Get the distance values from the read-out table and convert the strings into floats
min_MW_distance = [float(dist2[x, 0]) for x in range(len(dist2))]

#  get the velocity for a colormap
min_vel_mag = [trafo * np.sqrt(float(dist2[x, 9]) ** 2 + float(dist2[x, 10])
                               ** 2 + float(dist2[x, 11]) ** 2) for x in range(len(dist2))]

# extract galactic longitude
min_gal_l = [float(dist2[x, 4])
             for x in range(len(dist2)) if float(dist2[x, 1]) < 40]

# extract galactic latitude
min_gal_b = [float(dist2[x, 5])
             for x in range(len(dist2)) if float(dist2[x, 1]) < 40]

# now feed the coordinate lists into skycoord arrays
min_gal = SkyCoord(
    min_gal_l[:], min_gal_b[:], frame="galactic", unit="deg")


# do full length lists for trafo to Galactocentric system
# sun distance
min_sun_distance_trafo = [float(dist2[x, 1]) for x in range(len(dist2))]

# extract galactic longitude
min_gal_l_trafo = [float(dist2[x, 4]) for x in range(len(dist2))]

# extract galactic latitude
min_gal_b_trafo = [float(dist2[x, 5]) for x in range(len(dist2))]

# transformation to radial ccordinates around MW centre
min_x = [8.17827 - (min_sun_distance_trafo[i] * np.cos(min_gal_b_trafo[i] * np.pi / 180)
                    * np.cos(min_gal_l_trafo[i] * np.pi / 180)) for i in range(len(dist2))]

min_y = [min_sun_distance_trafo[i] *
         np.cos(min_gal_b_trafo[i] * np.pi / 180) * np.sin(min_gal_l_trafo[i] * np.pi / 180) for i in range(len(dist2))]

min_z = [min_sun_distance_trafo[i] *
         np.sin(min_gal_b_trafo[i] * np.pi / 180) for i in range(len(dist2))]

# min_x[0] = 1  # -370.7
# min_y[0] = 1  # 612.7
# min_z[0] = 1  # -283.1

min_r = [np.sqrt(min_x[i] ** 2 + min_y[i] **
                 2 + min_z[i] ** 2) for i in range(len(dist2))]

# and transform back to spherical for now Galactocentric coordinates
min_galcentre_l = [np.arctan2(min_y[i], min_x[i])
                   * (180 / np.pi) for i in range(len(dist2))]

min_galcentre_b = [np.arcsin(min_z[i] / min_r[i])
                   * (180 / np.pi) for i in range(len(dist2))]

# print(min_r[0])
# print(min_galcentre_l[0])
# print(min_galcentre_b[0])

# now feed the coordinate lists into skycoord arrays
min_galcentre = SkyCoord(
    min_galcentre_l[:], min_galcentre_b[:], frame="galactic", unit="deg")


#
'''PLOTTING:'''
#

# PRESENT TIME DATA SET-----------------------------------------------------------

# skymaps with distance colormap
plt.figure(1)
plt.subplot(111, projection='aitoff')
plt.grid(True)
# generate colourmap
cm = plt.cm.get_cmap('plasma')
# distance colormap
dist_cm = present_sun_distance
# plot the points as a scatter plot with distance cm
star_pos = plt.scatter(present_gal.l.wrap_at(
    '180d').radian, present_gal.b.radian, s=15, vmin=0, vmax=40, c=dist_cm, cmap=cm, alpha=0.7)

plt.colorbar(
    star_pos, label='$r_{Sun}$[kpc]', fraction=0.03, pad=0.03)
plt.title('Positions of simulation HVSs at $t_0$ \n Distance to Sun $r_{Sun}<40$ kpc \n Equal-mass scenario, Gal. coord.',
          loc='right', fontsize=10)
# plt.legend(loc='upper right')
plt.savefig("Present-star-dist-colourmap-same-mass.pdf")


'''
# plot with MW centre  in the middle
plt.figure(2)
plt.subplot(111, projection='aitoff')
plt.grid(True)
# generate colourmap
cm = plt.cm.get_cmap('plasma')
z = present_MW_distance
# plot the points as a scatter plot
star_pos = plt.scatter(present_galcentre.l.wrap_at('180d').radian,
                       present_galcentre.b.radian, s=15, vmin=0, vmax=50, c=z, cmap=cm, alpha=0.7)
plt.colorbar(
    star_pos, label='$r_{MW}$[kpc]', fraction=0.03, pad=0.04)
plt.title('Positions of simulation HVSs at $t_0$ \n Distance to MW centre $r_{MW}<50$ kpc \n Equal-mass scenario, Gal. coord.',
          loc='right', fontsize=10)
# plt.legend(loc='upper right')
# plt.savefig("Present-star-dist-colourmap-same-mass.pdf")
'''


# MINIMUM DISTANCE DATA SET-----------------------------------------------------------

# skymaps with colourmaps
# plot with SUn in the centre (Galactic coordinates)
plt.figure(3)
plt.subplot(111, projection='aitoff')
plt.grid(True)
# generate colourmap
cm = plt.cm.get_cmap('plasma')
z2 = min_sun_distance
# plot the points as a scatter plot
star_pos = plt.scatter(min_gal.l.wrap_at('180d').radian, min_gal.b.radian, s=15,\
            label='HVS positions', vmin=0, vmax=40, c=z2, cmap=cm, alpha=0.7)
plt.scatter(andromeda_gal.l.wrap_at('180d').radian, andromeda_gal.b.radian, s=100, \
            label='Andromeda position', marker='*', color='black')
plt.colorbar(star_pos, label='$r_{Sun}$[kpc]', fraction=0.03, pad=0.03)
plt.title(
    'Positions of simulation HVSs at minimum distance to MW centre \n Distance to Sun $r_{Sun}<40$ kpc \n Equal-mass scenario, Gal. coord.', loc='right', fontsize=10)
plt.legend(loc='lower right')
plt.savefig("Min-star-dist-colourmap-same-mass.pdf")
plt.show()

'''
# plot with MW centre in the middle
plt.figure(4)
plt.subplot(111, projection='aitoff')
plt.grid(True)
# generate colourmap
cm = plt.cm.get_cmap('plasma')
z2 = min_MW_distance
# plot the points as a scatter plot
star_pos = plt.scatter(min_galcentre.l.wrap_at('180d').radian,
                       min_galcentre.b.radian, s=15, vmin=0, vmax=50, c=z2, cmap=cm, alpha=0.7)
plt.colorbar(
    star_pos, label='$r_{MW}$ [kpc]', fraction=0.03, pad=0.04)
plt.title(
    'Positions of simulation HVSs at minimum distance to MW centre \n Distance to MW centre $r_{MW}<50$ kpc \n Equal-mass scenario, Gal. coord.', loc='right', fontsize=10)
#plt.legend(loc='upper right')
# plt.savefig("Min-star-dist-colourmap-same-mass.pdf")
plt.show()
'''
