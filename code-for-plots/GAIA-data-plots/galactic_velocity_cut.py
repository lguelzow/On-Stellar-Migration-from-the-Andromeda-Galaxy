import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

import numpy as np
import pandas as pd

import astropy
import healpy as hp
import sys
import scipy.stats
from astropy.coordinates import SkyCoord
import astropy.units as u
from mpl_toolkits.mplot3d import Axes3D


# CONSTANTS

# velocity for which to filter
filter_velocity = 0

# define constant for transformation from mas/yr to km/s
mas_trafo = 4.744213026

# define sun velocity vector
sun_vel = [11.1, 12.24 + 240, 7.25]
LSR = [0, 240, 0]

# limits of 3D graph axes
limits = 50


#
'''READ-OUT DATA FILES: '''
#

# readout process nr 1
# Read out the data of all minimum distances from the file
# comment in whichecer you need

data = pd.read_csv('/home/lguelzow/Nextcloud/MA/MA_Paper/GAIA/parallax-zero-point/Revision_GAIA_only_fast_test.csv')
# data = pd.read_csv('/home/lguelzow/Nextcloud/MA/MA_Paper/GAIA/190_25percent_error/error33_v_150-result.csv')
# data = pd.read_csv('/home/lguelzow/Nextcloud/MA/MA_Paper/GAIA/400kms-analysis/test2.csv')

print(data.columns)

print("Reading data into np.arrays...")

# read out data colums into np arrays
distance = np.array(data.parallax) ** (-1)

parallax = np.array(data.parallax)

parallax_zpt = np.array(data.parallax_zpt)

parallax_error = np.array(data.parallax_error)

# define corrected parallax
parallax_corr = parallax - parallax_zpt

# and parallax error criterion
parallax_criterion = parallax_error / parallax_corr

# print("Parallax \n", parallax[0:10])
# print("Parallax Zero Point \n", parallax_zpt[0:10])
# print("Corrected parallax \n", parallax_corr[0:10])
# print("Parallax error \n", parallax_error[0:10])
# print("Parallax criterion \n", parallax_criterion[0:10])


radial_velocity = np.array(data.radial_velocity)

radial_velocity_error = np.array(data.radial_velocity_error)

pm_ra = np.array(data.pmra)

pm_ra_error = np.array(data.pmra_error)

pm_dec = np.array(data.pmdec)

pm_dec_error = np.array(data.pmdec_error)

# print(pm_dec[0:100])

rig_asc = np.array(data.ra)

rig_asc_error = np.array(data.ra_error)

decl = np.array(data.dec)

decl_error = np.array(data.dec_error)

# define corrected parallax
parallax_corr = parallax - parallax_zpt

# print(pm_dec[0:10])

print("Finished! Preparing to filter data...")

# filtered lists for later
parallax_filter = []
parallax_zpt_filter = []
parallax_error_filter = []

rig_asc_filter = []
rig_asc_error_filter = []

decl_filter = []
decl_error_filter = []

radial_velocity_filter = []
radial_velocity_error_filter = []

pm_ra_filter = []
pm_ra_error_filter = []

pm_dec_filter = []
pm_dec_error_filter = []


# Define the velocity magnitude array
vel_mag = np.sqrt(radial_velocity ** 2 + ((pm_ra / parallax) * mas_trafo) ** 2 \
        + ((pm_dec / parallax) * mas_trafo) ** 2)




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
    phi_star = rig_asc[x] * (np.pi / 180) * (-1)
    theta_star = decl[x] * (np.pi / 180)  # * (-1)

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
gal_v = SkyCoord(vel_ra[:], vel_dec[:], frame="icrs",
                 unit="deg", equinox='J2015.5')


# convert to galactic coordinates via Skycoord
gal_v = gal_v.galactic

#
'''FILTER VELOCITIES IN GALACTO-CENTRIC FRAME: '''
#


# Calculation in galactocentric rest frame to filter out lower velocities:

# define list for kartesian gaia vector in galactic system and for filter of lower velocities
gaia_vel = [None] * (len(distance))

# list to store the changed velocity magnitude for filtering
vel_mag_criterium = [None] * (len(distance))
print("Starting filtering process...")

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
    

    # fill all filtered lists
    if vel_mag_criterium[i] > filter_velocity and distance[i] < limits and parallax_corr[i] > 0 and parallax_criterion[i] < 0.2:
        parallax_filter.append(parallax[i])
        parallax_zpt_filter.append(parallax_zpt[i])
        parallax_error_filter.append(parallax_error[i])

        rig_asc_filter.append(rig_asc[i])
        rig_asc_error_filter.append(rig_asc_error[i])

        decl_filter.append(decl[i])
        decl_error_filter.append(decl_error[i])

        radial_velocity_filter.append(radial_velocity[i])
        radial_velocity_error_filter.append(radial_velocity_error[i])

        pm_ra_filter.append(pm_ra[i])
        pm_ra_error_filter.append(pm_ra_error[i])

        pm_dec_filter.append(pm_dec[i])
        pm_dec_error_filter.append(pm_dec_error[i])


print("\n")
print(len(parallax_filter))


# produce a shorter file for long GAIA data files with only the fast stars
np.savetxt("Revision_GAIA_only_fast_testtest.csv", [])
f = open("Revision_GAIA_only_fast_testtest.csv", "w")

# write labels into first lines
f.write("Standin-for-sourceID" + "," + "parallax" + "," + "parallax_zpt" + "," + "ra" + "," + "dec" + "," + "pmra" + "," + "pmdec" + "," + "radial_velocity" + "," +
         "parallax_error" + "," + "ra_error" + "," + "dec_error" +
        "," + "pmra_error" + "," + "pmdec_error" + "," + "radial_velocity_error" + "\n")

# write into file the filtered lists
for i in range(len(parallax_filter)):
    f.write(str(0) + "," + str(parallax_filter[i]) + "," + str(parallax_zpt_filter[i]) + "," + str(rig_asc_filter[i]) + "," + str(decl_filter[i]) + "," +
            str(pm_ra_filter[i]) + "," + str(pm_dec_filter[i]) + "," + str(radial_velocity_filter[i]) + "," + str(parallax_error_filter[i]) +
              "," + str(rig_asc_error_filter[i]) + "," + str(decl_error_filter[i]) + "," +
            str(pm_ra_error_filter[i]) + "," + str(pm_dec_error_filter[i]) + "," + str(radial_velocity_error_filter[i]) + " " + "\n")
