import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
import pandas as pd

import astropy
import astropy.coordinates as coord
import healpy as hp
import sys
import scipy.stats
from astropy.coordinates import SkyCoord
import astropy.units as u
from mpl_toolkits.mplot3d import Axes3D


# function for healpy map
def cat2hpx(lon, lat, nside, radec=True):
    """
    Convert a catalogue to a HEALPix map of number counts per resolution
    element.

    Parameters
    ----------
    lon, lat : (ndarray, ndarray)
        Coordinates of the sources in degree. If radec=True, assume input is in the icrs
        coordinate system. Otherwise assume input is glon, glat

    nside : int
        HEALPix nside of the target map

    radec : bool
        Switch between R.A./Dec and glon/glat as input coordinate system.

    Return
    ------
    hpx_map : ndarray
        HEALPix map of the catalogue number counts in Galactic coordinates

    """

    npix = hp.nside2npix(nside)

    if radec:
        eq = SkyCoord(lon, lat, frame="icrs", unit='deg')
        l, b = eq.galactic.l.value, eq.galactic.b.value
    else:
        l, b = lon, lat

    # conver to theta, phi
    theta = np.radians(90. - b)
    phi = np.radians(l)

    # convert to HEALPix indices
    indices = hp.ang2pix(nside, theta, phi)

    idx, counts = np.unique(indices, return_counts=True)

    # fill the fullsky map
    hpx_map = np.zeros(npix, dtype=int)
    hpx_map[idx] = counts

    return hpx_map


# CONSTANTS

# velocity for which to filter
filter_velocity = 500

# define constant for transformation from mas/yr to km/s
mas_trafo = 4.744213026

# transformation between km/s and kpc/Myr
trafo = 1.023 / 1000  # 1 kpc/Myr = (1000 / 1.023) km/s

# define sun velocity vector
sun_vel = [11.1, 12.24 + 240, 7.25]
LSR = [0, 240, 0]

# x coordinate of Galactic centre in Galactic coordiates
x_GC = 8.178

# scaling factor for size of velocity vectors in 3D plot
scale = 3

# limits of 3D graph axes
limits = 40

limits_plot = 25
# 7

# Andromeda position in galactic coordinates
andromeda_pos = [121.172, -21.573]

andromeda_gal = SkyCoord(andromeda_pos[0], andromeda_pos[1], frame="galactic", unit="deg")


#
'''READ-OUT DATA FILES: '''
#

# readout process nr 1
# Read out the data of all minimum distances from the file
# comment in whichecer you need

# data = pd.read_csv('/home/lguelzow/Nextcloud/MA/MA_Paper/GAIA/190_25percent_error/Revision_GAIA_only_fast_33.csv')
data = pd.read_csv('/home/lguelzow/Nextcloud/MA/MA_Paper/GAIA/Revision_GAIA_only_fast_test.csv')
# data = pd.read_csv('/home/lguelzow/Nextcloud/MA/MA_Paper/GAIA/145kms-new-filter/Paper_GAIA_only_fast_500.csv')
# data = pd.read_csv('/home/lguelzow/Nextcloud/MA/MA_Paper/GAIA/190kms-analysis/Paper_GAIA_only_fast_500_2.csv')

print(data.columns)

print("Reading data into np.arrays...")

# read out data colums into np arrays
# distance = np.array(data.parallax) ** (-1)

parallax = np.array(data.parallax)

parallax_zpt = np.array(data.parallax_zpt)

# define corrected parallax
parallax_corr = parallax - parallax_zpt

distance = parallax_corr ** (-1)

radial_velocity = np.array(data.radial_velocity)

pm_ra = np.array(data.pmra)

pm_dec = np.array(data.pmdec)

# print(pm_dec[0:100])

right_ascension = np.array(data.ra)

declination = np.array(data.dec)

# print(pm_dec[0:10])

print("Finished! Preparing to filter data...")

# filtered lists for later
distance_filter = []
parallax_filter =[]


# Define the velocity magnitude array
vel_mag = np.sqrt(radial_velocity ** 2 + ((pm_ra / parallax) * mas_trafo) ** 2 \
        + ((pm_dec / parallax) * mas_trafo) ** 2)

# as well as filtered list for later
vel_mag_filter = []


# Define the right ascension and declination arrays
# as well as filtered list for later
rig_asc = np.array(data.ra)
decl = np.array(data.dec)

# as well as filtered lists for later
rig_asc_filter = []
decl_filter = []


# converting the RA_Dec coordinates into Skycoord for plotting and transformations
radec = SkyCoord(rig_asc[:], decl[:], frame="icrs",
               unit="deg", equinox='J2015.5')

# transformation to galactic coordinates
gal = radec.galactic


# get kartesian positions of stars (galactic system)
# list for position vectors
# as well as filtered list for later
kart_posi = np.array([[distance[x] * np.cos(gal.b.deg[x] * (np.pi / 180)) * np.cos(gal.l.deg[x] * (np.pi / 180)), \
                       distance[x] * np.cos(gal.b.deg[x] * (np.pi / 180)) * np.sin(gal.l.deg[x] * (np.pi / 180)),\
                       distance[x] * np.sin(gal.b.deg[x] * (np.pi / 180))] for x in range(len(distance))])

# kart_posi = [None] * (len(distance))
kart_posi_filter = []


# Define the star velocity array while excluding the 1st line
# and arrays for the components of their spherical coordinates
vel_ra = [None] * (len(distance))
vel_dec = [None] * (len(distance))

# Get the radial velocity and proper motion values values from the read-out table and convert the strings into floats
# while also converting them into the real velocities
vel_vec = np.array([[radial_velocity[x], (pm_ra[x] / parallax[x]) * mas_trafo,
                  (pm_dec[x] / parallax[x]) * mas_trafo] for x in range(len(distance))])


for x in range(len(distance)):

    # record progress of calculation
    if x % 100 == 0:
        sys.stdout.write("\rProgress-Trafo:" + str(float(x) / len(distance)))
        sys.stdout.flush()

    # vel_vec[x] = np.array([radial_velocity[x], (pm_ra[x] / parallax[x]) * mas_trafo,
    #               (pm_dec[x] / parallax[x]) * mas_trafo])

    # if conditions for testing purposes
    # if conditions to delete transverse velocity and to make all radial velocity positive
    # if abs(pm_ra[x]) < 0:
    #     vel_vec[x][1] = 0
    # if abs(pm_dec[x]) < 0:
    #     vel_vec[x][2] = 0
    # if float(gaia[x + 1, 6]) < 0:
    #    vel_vec[x][0] = vel_vec[x][0] * (-1)
    # print(vel_vec[51])

    # Define rotation matrix to rotate velocity direction so that it's right in the global coordinate system and not the star's
    # First get the angles to rotate around from the position of the star in radians
    phi_star = right_ascension[x] * (np.pi / 180) * (-1)
    theta_star = declination[x] * (np.pi / 180)  # * (-1)

    # rotation matrix around z-axis with RA angle of star
    Rot_z = np.array([[np.cos(phi_star), np.sin(phi_star), 0],
                      [np.sin(phi_star) * (-1), np.cos(phi_star), 0], [0, 0, 1]])

    # rotation around y-axis with DEC angle of star
    Rot_y = np.array([[np.cos(theta_star), 0, np.sin(theta_star) * (-1)], [0, 1, 0], [
                     np.sin(theta_star), 0, np.cos(theta_star)]])

    # perform both rotations
    vel_vec[x] = np.dot(Rot_y, vel_vec[x])
    vel_vec[x] = np.dot(Rot_z, vel_vec[x])

    # transform to spherical coordinates = [RA, DEC]
    vel_ra[x] = np.arctan2(vel_vec[x][1], vel_vec[x][0]) * (180 / np.pi)
    vel_dec[x] = np.arcsin(vel_vec[x][2] / (vel_vec[x][0] ** 2 +
                                            vel_vec[x][1] ** 2 + vel_vec[x][2] ** 2) ** 0.5) * (180 / np.pi)


#
''' ALTERNATIVE METHOD: '''
#
'''
vel_ra_alt = [None] * (len(gaia) - 1)
vel_dec_alt = [None] * (len(gaia) - 1)

# Define the array for the kartesian velocity vector
vel_kart = [None] * (len(gaia) - 1)

# Get the right ascension values from the read-out table and convert the strings into floats
for x in range(len(gaia) - 1):
    # establish as vector
    vel_kart[x] = [None, None, None]

    # readout angles and convert to radian
    RA = float(gaia[x + 1, 2]) * (np.pi / 180)#  * (-1)
    DEC = (90 - float(gaia[x + 1, 3])) * (np.pi / 180)

    # x-component: dx/dt
    vel_kart[x][0] = float(gaia[x + 1, 6]) * np.sin(DEC) * np.cos(RA) - (float(gaia[x + 1, 5]) / float(gaia[x + 1, 1])) * mas_trafo * \
        np.cos(DEC) * np.cos(RA) - (float(gaia[x + 1, 4]) / float(
            gaia[x + 1, 1])) * mas_trafo * np.sin(DEC) * np.sin(RA) # / np.cos(float(gaia[x + 1, 3]) * (np.pi / 180))

    # y-component: dy/dt
    vel_kart[x][1] = float(gaia[x + 1, 6]) * np.sin(DEC) * np.sin(RA) - (float(gaia[x + 1, 5]) / float(gaia[x + 1, 1])) * mas_trafo * \
        np.cos(DEC) * np.sin(RA) + (float(gaia[x + 1, 4]) / float(
            gaia[x + 1, 1])) * mas_trafo * np.sin(DEC) * np.cos(RA) # / np.cos(float(gaia[x + 1, 3]) * (np.pi / 180))

    # z-component: dz/dt
    vel_kart[x][2] = float(gaia[x + 1, 6]) * np.cos(DEC) + (
        float(gaia[x + 1, 5]) / float(gaia[x + 1, 1])) * mas_trafo * np.sin(DEC)

    # transform to spherical coordinates = [RA, DEC]
    vel_ra_alt[x] = np.arctan2(vel_kart[x][1], vel_kart[x][0]) * (180 / np.pi)
    vel_dec_alt[x] = np.arcsin(vel_kart[x][2] / (vel_kart[x][0] ** 2 +
                                                 vel_kart[x][1] ** 2 + vel_kart[x][2] ** 2) ** 0.5) * (180 / np.pi)


# test for left vs right-handed
# Get the right ascension values from the read-out table and convert the strings into floats
for x in range(len(gaia) - 1):
    # establish as vector
    vel_kart[x] = [None, None, None]

    # readout angles and convert to radian
    RA = float(gaia[x + 1, 2]) * (np.pi / 180) # * (-1)
    DEC = float(gaia[x + 1, 3]) * (np.pi / 180)

    # x-component: dx/dt
    vel_kart[x][0] = float(gaia[x + 1, 6]) * np.cos(DEC) * np.cos(RA) - (float(gaia[x + 1, 5]) / float(gaia[x + 1, 1])) * mas_trafo * \
        np.sin(DEC) * np.cos(RA) - (float(gaia[x + 1, 4]) / float(
            gaia[x + 1, 1])) * mas_trafo * np.cos(DEC) * np.sin(RA) # / np.cos(float(gaia[x + 1, 3]) * (np.pi / 180))

    # y-component: dy/dt
    vel_kart[x][1] = float(gaia[x + 1, 6]) * np.cos(DEC) * np.sin(RA) - (float(gaia[x + 1, 5]) / float(gaia[x + 1, 1])) * mas_trafo * \
        np.sin(DEC) * np.sin(RA) + (float(gaia[x + 1, 4]) / float(
            gaia[x + 1, 1])) * mas_trafo * np.cos(DEC) * np.cos(RA) # / np.cos(float(gaia[x + 1, 3]) * (np.pi / 180))

    # z-component: dz/dt
    vel_kart[x][2] = float(gaia[x + 1, 6]) * np.sin(DEC) + (
        float(gaia[x + 1, 5]) / float(gaia[x + 1, 1])) * mas_trafo * np.cos(DEC)

    # transform to spherical coordinates = [RA, DEC]
    vel_ra_alt[x] = np.arctan2(vel_kart[x][1], vel_kart[x][0]) * (180 / np.pi)
    vel_dec_alt[x] = np.arcsin(vel_kart[x][2] / (vel_kart[x][0] ** 2 +
                                                 vel_kart[x][1] ** 2 + vel_kart[x][2] ** 2) ** 0.5) * (180 / np.pi)

# print(vel_kart[202])
# print(vel_ra_alt[202])
# print(vel_dec_alt[202])


# import coordinates in astropy funktion "Skycoord"
gal_v = SkyCoord(vel_ra_alt[:], vel_dec_alt[:], frame="icrs", unit="deg")

'''
# import coordinates in astropy funktion "Skycoord"
gal_v = SkyCoord(vel_ra[:] * u.degree, vel_dec[:] * u.degree, frame="icrs",
                 equinox='J2015.5', distance=distance * u.kpc)


# convert to galactic coordinates via Skycoord
gal_v = gal_v.galactic

#
'''FILTER VELOCITIES IN GALACTO-CENTRIC FRAME: '''
#


# Calculation in galactocentric rest frame to filter out lower velocities:

# define list for kartesian gaia vector in galactic system and for filter of lower velocities
gaia_vel = [None] * (len(distance))
gaia_vel_filter = []

# filtered galactic coordinates list
vel_ra_filter = []
vel_dec_filter = []

# lists for filtered output data set
pmra_filter = []
pmdec_filter = []
v_rad_filter = []

# galactic_vel_l = [None] * (len(distance))
# galactic_vel_b = [None] * (len(distance))

gal_rest_vel_filter = []

# list to store the changed velocity magnitude for filtering
vel_mag_criterium = [None] * (len(distance))

for i in range(len(distance)):

    # record progress of calculation
    if i % 100 == 0:
        sys.stdout.write("\rProgress-Filter:" +
                         str(float(i) / (len(distance))))
        sys.stdout.flush()

    # construct velocity vector from gal. coord. and subtract sun velocity
    gaia_vel[i] = [vel_mag[i] * np.cos(gal_v.b.radian[i]) * np.cos(gal_v.l.radian[i]), vel_mag[i]
                   * np.cos(gal_v.b.radian[i]) * np.sin(gal_v.l.radian[i]), vel_mag[i] * np.sin(gal_v.b.radian[i])]

    # add sun velocity vector to recalculate velocity magnitude in galactocentric frame for filtering
    vel_mag_criterium[i] = ((gaia_vel[i][0] + sun_vel[0]) ** 2 + (
        gaia_vel[i][1] + sun_vel[1]) ** 2 + (gaia_vel[i][2] + sun_vel[2]) ** 2) ** 0.5
    
    sun_distance = np.sqrt(kart_posi[i][0] ** 2 + kart_posi[i][1] ** 2 + kart_posi[i][2] ** 2)

    # fill all filtered lists
    if vel_mag_criterium[i] > filter_velocity and sun_distance < limits: # distance[i] < limits:
        parallax_filter.append(parallax[i])
        distance_filter.append(distance[i])
        rig_asc_filter.append(rig_asc[i])
        decl_filter.append(decl[i])
        gaia_vel_filter.append(gaia_vel[i])
        vel_mag_filter.append(vel_mag[i])
        kart_posi_filter.append(kart_posi[i])
        vel_ra_filter.append(gal_v.l.deg[i])
        vel_dec_filter.append(gal_v.b.deg[i])
        gal_rest_vel_filter.append(vel_mag_criterium[i])

        # for output data set to reproduce the format of the original data
        pmra_filter.append(pm_ra[i])
        pmdec_filter.append(pm_dec[i])
        v_rad_filter.append(radial_velocity[i])

# transform to spherical coordinates = [gal_l, gal_b]
galactic_vel_l = np.array([np.arctan2(gaia_vel[i][1] + sun_vel[1], gaia_vel[i][0] + sun_vel[0]) * (180 / np.pi) for i in range(len(distance))])
galactic_vel_b = np.array([np.arcsin((gaia_vel[i][2] + sun_vel[2]) / vel_mag_criterium[i]) * (180 / np.pi) for i in range(len(distance))])

# heliocentric velocity
helio_vel_l = np.array([np.arctan2(gaia_vel[i][1], gaia_vel[i][0]) * (180 / np.pi) for i in range(len(distance))])
helio_vel_b = np.array([np.arcsin((gaia_vel[i][2]) / vel_mag[i]) * (180 / np.pi) for i in range(len(distance))])

# import velocity vector in Galactocentric rest frame into Skycoord
gal_v_rest_frame = SkyCoord(galactic_vel_l[:], galactic_vel_b[:], frame="galactic", unit="deg")

kart_posi_filter = np.array(kart_posi_filter)

# radec_v_rest_frame = gal_v_rest_frame.icrs
# print(radec_v_rest_frame)

print("\n")
print("Amount of analysed stars: ", len(parallax_filter))

# list for kartesian vector plot and combined list
kart_vel = [None] * (len(gaia_vel_filter))
kart_trajectory = np.array([None] * (len(gaia_vel_filter)))

kart_trajectory_gal_rest = np.array([[kart_posi_filter[i][0], kart_posi_filter[i][1], kart_posi_filter[i][2], \
                                     (gaia_vel[i][0] + sun_vel[0])  * trafo * scale, (gaia_vel[i][1] + sun_vel[1])  * trafo * scale, (gaia_vel[i][2] + sun_vel[2]) * trafo * scale] 
                                     for i in range(len(gal_rest_vel_filter))])
print(kart_trajectory_gal_rest)

# lists for velocities in Toomre diagram of velocities
toomre_vel_perpendicular = [None] * (len(gaia_vel_filter))
toomre_vel_disc = [None] * (len(gaia_vel_filter))

# lists for velocity component histograms
x_hist_vel = [None] * (len(gaia_vel_filter))
y_hist_vel = [None] * (len(gaia_vel_filter))
z_hist_vel = [None] * (len(gaia_vel_filter))

# average velocity vector direction
av_vel = [0, 0, 0]

# filtered galactic coordinates list
vel_ra_filter2 = [None] * (len(gaia_vel_filter))
vel_dec_filter2 = [None] * (len(gaia_vel_filter))

# transform filtered data to galactic coordinates again
# loop for Toomre diagram and component histograms with filtered data
for x in range(len(gaia_vel_filter)):
    # code for Toomre diagram and component histograms
    # fill list of kartesian vel. vectors in km/s, already in galacto-centric rest frame
    kart_vel[x] = gaia_vel_filter[x]

    # add up all velocities
    av_vel[0] += gaia_vel_filter[x][0]
    av_vel[1] += gaia_vel_filter[x][1]
    av_vel[2] += gaia_vel_filter[x][2]

    # fill lists for Toomre diagram
    toomre_vel_perpendicular[x] = (
        kart_vel[x][0] ** 2 + kart_vel[x][2] ** 2) ** 0.5
    toomre_vel_disc[x] = kart_vel[x][1]

    # fill lists for component histograms
    x_hist_vel[x] = kart_vel[x][0]
    y_hist_vel[x] = kart_vel[x][1]
    z_hist_vel[x] = kart_vel[x][2]

    # transform to kpc/Myr for plot and enlarge the vectors to make them more visible
    kart_vel[x][0] = kart_vel[x][0] * trafo * scale
    kart_vel[x][1] = kart_vel[x][1] * trafo * scale
    kart_vel[x][2] = kart_vel[x][2] * trafo * scale

    # and combine both kartesian lists
    kart_trajectory[x] = [kart_posi_filter[x][0], kart_posi_filter[x][1], kart_posi_filter[x][2],
                          kart_vel[x][0], kart_vel[x][1], kart_vel[x][2]]

    # transform back to gal. coord.
    vel_ra_filter2[x] = np.array(np.arctan2(
        gaia_vel_filter[x][1], gaia_vel_filter[x][0]) * (180 / np.pi))
    vel_dec_filter2[x] = np.array(np.arcsin(gaia_vel_filter[x][2] / (gaia_vel_filter[x][0] ** 2 +
                                                            gaia_vel_filter[x][1] ** 2 + gaia_vel_filter[x][2] ** 2) ** 0.5) * (180 / np.pi))

# and divide by number of them to get average velocity vector
av_vel[0] = av_vel[0] / len(gaia_vel_filter)
av_vel[1] = av_vel[1] / len(gaia_vel_filter)
av_vel[2] = av_vel[2] / len(gaia_vel_filter)

# transform to Galactic coordinates
av_vel_gal_l = np.arctan2(av_vel[1], av_vel[0]) * (180 / np.pi)
av_vel_gal_b = np.arcsin(
    av_vel[2] / (av_vel[0] ** 2 + av_vel[1] ** 2 + av_vel[2] ** 2) ** 0.5) * (180 / np.pi)

# put into Skycoord
av_vel_gal = SkyCoord(av_vel_gal_l, av_vel_gal_b, frame='galactic', unit='deg')

print(av_vel_gal)
# print(vel_dec_filter[0:10])



#
# PLOTTING:
#

# generate fit for velocity magnitude histogram
# mu, sigma = scipy.stats.expon.fit(vel_mag)
# best_fit_line = scipy.stats.expon.pdf(100, mu, sigma)
# print(sigma)

# heliocentric velocity of Sgr A*
sgr_disc = -sun_vel[1]
sgr_perp = (sun_vel[0] ** 2 + sun_vel[2] ** 2) ** 0.5

# Toomre diagram of star velocities
plt.figure(1)
toomre = plt.gca()
toomre.scatter(toomre_vel_disc, toomre_vel_perpendicular,
               2, label='Gaia HVS ($d<50$kpc)')
toomre.scatter(sgr_disc, sgr_perp, 15, 'lime', label='Sgr A*')
toomre.scatter(0, 0, 15, 'crimson', label='Sun')
# add circles of equal velocity to plot
circle_sun = plt.Circle((-252, 13.26), 252, color='crimson',
                        fill=False)
circle_min = plt.Circle((sgr_disc, sgr_perp), 600, color='black', fill=False)
circle_mid = plt.Circle((sgr_disc, sgr_perp), 700, color='black', fill=False)
circle_high = plt.Circle((sgr_disc, sgr_perp), 800,
                         color='black', fill=False)
circle_max = plt.Circle((sgr_disc, sgr_perp), 1000,
                        color='black', fill=False)
circle_out = plt.Circle((sgr_disc, sgr_perp), 1500, color='black', fill=False)
toomre.add_artist(circle_sun)
toomre.add_artist(circle_min)
toomre.add_artist(circle_mid)
toomre.add_artist(circle_high)
toomre.add_artist(circle_max)
# toomre.add_artist(circle_out)
# formatting
toomre.set_xlabel('$v_y$ [km/s]')
toomre.set_ylabel(r'$\sqrt{(v^2_x + v^2_z)}$ [km/s]')
toomre.set_xlim(-1200, 1000)
toomre.set_ylim(0, 1000 + sgr_perp)
toomre.set_aspect('equal')
toomre.legend(loc='upper right')
plt.savefig("GAIA-Vel_Toomre_diagram.pdf")


random_selection = np.full(len(kart_trajectory), False)
random_selection[:100] = True

np.random.shuffle(random_selection)
print(random_selection)

# kartesian plot of star positions and velocities
# for 13kpc  ([300:450])
# X, Y, Z, U, V, W = zip(*kart_trajectory[random_selection])
X, Y, Z, U, V, W = zip(*kart_trajectory_gal_rest[random_selection])
# for 5kpc
# X, Y, Z, U, V, W = zip(*kart_trajectory[100:150])

# repeat for radec vectors
# X1, Y1, Z1, U1, V1, W1 = zip(*kart_trajectory_radec[0:150])

# plot position of Sun and Milky Way
x_point_sun = 0
x_point_sgr = x_GC # 8.178
y_points = 0
z_points = 0


# galactic coordinates 3D plot
plt.figure(2)
ax = plt.subplot(111, projection='3d')
ax.quiver(X, Y, Z, U, V, W, label='Stars with $v_{GC}>400$km/s', alpha=1)

# plot positions of Sun and Sgr A*
ax.scatter(x_point_sun, y_points, z_points,
           c='crimson', s=40, label='Sun')
ax.scatter(x_point_sgr, y_points, z_points,
           c='green', s=40, label='Sgr A*')

ax.set_xlim([-limits_plot, limits_plot])
ax.set_ylim([-limits_plot, limits_plot])
ax.set_zlim([-limits_plot, limits_plot])
ax.set_xlabel('x-axis [kpc]')
ax.set_ylabel('y-axis [kpc]')
ax.set_zlabel('z-axis [kpc]')
ax.legend(loc='upper right')
plt.savefig("GAIA-3D.pdf")


fig = matplotlib.pyplot.figure(12)
ax = fig.add_subplot(111, projection='3d')
# ax= plt.subplots(111, projection='3d')
'''
ax.quiver(X1, Y1, Z1, U1, V1, W1, color='black', label=r'10.0 $\leq t$ < 10.5')

ax.quiver(X2, Y2, Z2, U2, V2, W2, color='navy', label=r'10.5 $\leq t$ < 11.0')
ax.quiver(X3, Y3, Z3, U3, V3, W3, color='dodgerblue',
          label=r'11.0 $\leq t$ < 11.5')
ax.quiver(X4, Y4, Z4, U4, V4, W4, color='turquoise',
          label=r'11.5 $\leq t$ < 12.0')
'''
# generate colourmap
cm = matplotlib.cm.coolwarm
norm = matplotlib.colors.LogNorm(vmin=400, vmax=1000)
c = gal_rest_vel_filter
sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
sm.set_array([])
# print(cm(c))
# print(c[0:30])
# 3D plot with colormap
q = ax.quiver(X, Y, Z, U, V, W, color=cm(norm(c)),
              label='Gaia HVSs with $v_{GC}>400$km/s', alpha=1, normalize=True, length=2.5)
# q.set_array(a)
plt.colorbar(sm, label="Velocity in Galactic rest frame [km/s]")

# plot positions of Sun and Sgr A*
ax.scatter(x_point_sun, y_points, z_points, c='black', s=40, label='Sun')
ax.scatter(x_point_sgr, y_points, z_points,
           c='lime', s=40, label='Sgr A*')


ax.set_xlim([-limits_plot, limits_plot])
ax.set_ylim([-limits_plot, limits_plot])
ax.set_zlim([-limits_plot, limits_plot])
ax.set_xlabel('$x$-axis [kpc]')
ax.set_ylabel('$y$-axis [kpc]')
ax.set_zlabel('$z$-axis [kpc]')
ax.legend(loc='upper right', fontsize=9)
plt.show()
plt.close()


'''
# radec 3D plot
plt.figure(18)
ax = plt.subplot(111, projection='3d')
ax.quiver(X1, Y1, Z1, U1, V1, W1)
ax.set_xlim([-limits, limits])
ax.set_ylim([-limits, limits])
ax.set_zlim([-limits, limits])
ax.set_xlabel('x-axis [kpc]')
ax.set_ylabel('y-axis [kpc]')
ax.set_zlabel('z-axis [kpc]')
# plt.savefig("GAIA-3D_radec.pdf")
'''
plt.show()
plt.close()


# HISTOGRAMS


# plot velocity magnitudes as a histogram
plt.figure(3)
# plt.hist(test_plot, 30)
# plt.hist(vel_mag_criterium, bins=25,  weights=np.ones(len(vel_mag_criterium)) / len(vel_mag_criterium))
plt.hist(vel_mag, bins=20) # , weights=np.ones(len(vel_mag_filter)) / len(vel_mag_filter))
plt.xlabel('Heliocentric velocity $v$ [km/s]')
plt.ylabel('# of HVSs')
plt.yscale('log')
# plt.xlim(400, 1700)
# plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.savefig('GAIA-Vel_mag_histogram.pdf')


plt.figure(33)
# plt.hist(test_plot, 30)
plt.hist(vel_mag_criterium, bins=25) # ,  weights=np.ones(len(vel_mag_criterium)) / len(vel_mag_criterium))
# plt.hist(vel_mag, bins=20) # , weights=np.ones(len(vel_mag_filter)) / len(vel_mag_filter))
plt.xlabel('Galactocentric velocity $v$ [km/s]')
plt.ylabel('# of HVSs')
plt.yscale('log')
# plt.xlim(400, 1700)
# plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.savefig('GAIA-Vel_mag_histogram_Gal.pdf')


# plot histograms for velocity components
# v_x radial component
plt.figure(4)
plt.hist(x_hist_vel, 20)# , weights=np.ones(len(x_hist_vel)) / len(x_hist_vel))
plt.xlabel('$v_x$ (radial) [km/s]')
plt.ylabel('Fraction of HVSs')
# plt.yscale("log")
# plt.xlim(500, 1000)
# plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.savefig("GAIA-Vel_x-comp_hist.pdf")

# v_y disc component
plt.figure(5)
plt.hist(y_hist_vel, 20) # , weights=np.ones(len(y_hist_vel)) / len(y_hist_vel))
plt.xlabel('$v_y$ (rotation) [km/s]')
plt.ylabel('Fraction of HVS')
# plt.yscale("log")
# plt.xlim(500, 1000)
# plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.savefig("GAIA-Vel_y-comp_hist.pdf")

# v_z pole component
plt.figure(6)
plt.hist(z_hist_vel, 20) # , weights=np.ones(len(z_hist_vel)) / len(z_hist_vel))
plt.xlabel('$v_z$ (polar) [km/s]')
plt.ylabel('Fraction of HVSs')
# plt.yscale("log")
# plt.xlim(500, 1000)
# plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.savefig("GAIA-Vel_z-comp_hist.pdf")


plt.figure(7)
# plt.hist(test_plot, 30)
plt.hist(distance, bins=20) # ,  weights=np.ones(len(distance)) / len(distance))
plt.xlabel('Distance to the Sun [kpc]')
plt.ylabel('Fraction of HVSs')
# plt.yscale('log')
# plt.xlim(0, 33)
# plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.savefig('GAIA-Minimum-distance_histogram.pdf')


# SKYMAPS


# plot distribution of results in the sky as a skymap
plt.figure(9)
ax1 = plt.subplot(111, projection='aitoff')
ax1.grid(True)

# generate colourmap
cm = plt.cm.get_cmap('plasma')
z = distance

star_pos = ax1.scatter(gal.l.wrap_at('180d').radian,
                       gal.b.radian, s=2, c=z, cmap=cm, alpha=0.7)
plt.colorbar(star_pos, label='Distance to Sun [kpc]',
             fraction=0.03, pad=0.04)
ax1.set_title('Positions of HVSs from GAIA DR3 measurements \n ($d < 50$ kpc, $v_{LSR} > 600$ km/s, Gal. coord.)',
              loc='right', fontsize=10)
plt.savefig("GAIA-Star_distribution_skymap_heatmap.pdf")


# Projection plot of star positions within the Milky Way
plt.figure(91)
plt.rcParams['font.size'] = 12
plt.title('Positions of Gaia stars with $v_\mathrm{GC}>500\,\mathrm{kms}^{-1}$\n Distance to Sun $r_\mathrm{Sun}<40\,\mathrm{kpc}$\n Projected onto the Galactic disc')
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.03)
# generate colourmap
# colormap
cm = matplotlib.cm.plasma
'''norm = matplotlib.colors.LogNorm(vmin=500, vmax=1000)
c = gal_rest_vel_filter
sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
sm.set_array([])'''
c = kart_posi_filter[:, 2]
# main scatter plot
proj_stars = ax.scatter(kart_posi_filter[:, 1], kart_posi_filter[:, 0] - x_GC,
               s=3, c=c, cmap=cm, vmin=-35, vmax=35) # , label=r'Gaia stars ($v_{GC}>500\,\mathrm{km/s})$')
ax.scatter(0, 0, s=50, color='black', label='Galactic Centre', marker='*')
ax.scatter(0,  -x_GC, s=40, color='black', label='Sun', marker='^')
# add circles of equal velocity to plot
circle_sun1 = plt.Circle((0, 0), 10, color='black', fill=False, linestyle='--')
circle_inner = plt.Circle((0, 0), 25, color='black', fill=False, linestyle='--')
# circle_mid = plt.Circle((0, 0), 30, color='black', fill=False, linestyle='--')
circle_outer = plt.Circle((0, 0), 50, color='black', fill=False, linestyle='--')
ax.add_artist(circle_sun1)
ax.add_artist(circle_inner)
# ax.add_artist(circle_mid)
ax.add_artist(circle_outer)
# add labels to circles
def text(x, y, text):
    ax.text(x, y, text,
            ha='left', va='center', color='black', fontsize=8.5, weight='bold')
text(-40.5, 27, '50 kpc')
# text(-24.5, 16.333333, '30 kpc')
text(-20, 13.333333, '25 kpc')
text(-7.5, 5, '10 kpc')

# add colorbar
plt.colorbar(proj_stars, label='Elevation above Galactic plane [kpc]', cax=cax) # fraction=0.03, pad=0.04) # , ticks=[500, 600, 700, 800, 900, 1000])
# formatting
ax.set_xlabel('Distance [kpc]')
ax.set_ylabel('Distance [kpc]')
ax.set_xlim(-51, 51)
ax.set_ylim(-34, 34)
ax.set_aspect('equal')
ax.legend(loc='upper right')
plt.savefig("GAIA-Star-pos-projection.pdf")



'''# same plot but as healpy map
NSide = 2**12
NPix= hp.nside2npix(NSide)
map=np.zeros(NPix,float)
hp.mollview(map,xsize=NSide)'''

hpx_map = cat2hpx(rig_asc_filter[:], decl_filter[:], nside=2 ** 5, radec=True)
hp.mollview(hpx_map+1, norm='log', unit='# of stars in bin (log. scale)', \
     title='Spatial distribution of stars from GAIA DR3 in Gal. coordinates\n  with $v_\mathrm{GC}>500\,\mathrm{kms}^{-1}$ and $r_\mathrm{Sun}<40\,\mathrm{kpc}$', cmap='plasma')
plt.rcParams['font.size'] = 12
plt.savefig("GAIA-Mollview.pdf")
plt.show()
plt.close()


# velocity skymap
# use galactic_vel_() for velocity in Galactocentric restframe
# use vel_ra_filter2 / vel_dec_filter2 for velocity in heliocentric restframe

hpx_map = cat2hpx(galactic_vel_l[:], galactic_vel_b[:], nside=2 ** 5, radec=False)
# hpx_map = cat2hpx(helio_vel_l[:], helio_vel_b[:], nside=2 ** 5, radec=False)
hp.mollview(hpx_map+1, flip='geo', norm='log', unit='# of stars in bin (log. scale)', \
     title='Velocity direction distribution of stars from GAIA DR3 in Gal. coordinates\n with $v_\mathrm{GC}>500\,\mathrm{kms}^{-1}$ and $r_\mathrm{GC}<50\,\mathrm{kpc}$', cmap='plasma')
plt.rcParams['font.size'] = 12
plt.savefig("GAIA-Mollview_vel.pdf")

# transform to spherical coordinates = [gal_l, gal_b]
galactic_vel_l = np.array([np.arctan2(gaia_vel[i][1] + sun_vel[1], gaia_vel[i][0] + sun_vel[0]) * (180 / np.pi) for i in range(len(distance))])
galactic_vel_b = np.array([np.arcsin((gaia_vel[i][2] + sun_vel[2]) / vel_mag_criterium[i]) * (180 / np.pi) for i in range(len(distance))])

galactic_rest_v = SkyCoord(galactic_vel_l[:] * u.degree, galactic_vel_b[:] * u.degree, frame="galactic")

# gal_v.transform_to(coord.Galactocentric(galcen_distance=0))
# gal_v = gal_v.galactic

# plot distribution of velocity directions at minimum distance as a skymap
plt.figure(11)
ax2 = plt.subplot(111, projection='aitoff', label='velocity')
ax2.grid(True)
# colormap
cm = matplotlib.cm.plasma
norm = matplotlib.colors.LogNorm(vmin=500, vmax=2500)
c = vel_mag_criterium
sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
sm.set_array([])

star_vel = ax2.scatter(gal_v.l.wrap_at('180d').radian,
                       gal_v.b.radian, s=5, color=cm(norm(c)), alpha=0.7, label=None)
# ax2.scatter(gal_v_rest_frame.l.wrap_at('180d').radian, gal_v_rest_frame.b.wrap_at('180d').radian, s=2, color=cm(norm(c)), alpha=0.7, label=None)
plt.colorbar(sm, label='Galactocentric Velocity [km/s]', fraction=0.03, pad=0.04)
ax2.legend(loc='upper right')
ax2.set_title(
    'Velocity directions of HVSs from GAIA DR3 measurements \n ($d < 50$ kpc, $v_{LSR} > 600$ km/s, Gal. coord.)', loc='right', fontsize=10)
plt.savefig("GAIA-Vel_distribution_skymap_heatmap.pdf")
plt.show()
plt.close()
# l.wrap_at('180d').