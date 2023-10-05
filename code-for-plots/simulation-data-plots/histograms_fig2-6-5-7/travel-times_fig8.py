# from matplotlib.lines import _LineStyle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib.ticker import PercentFormatter
from mpl_toolkits import mplot3d
import matplotlib.colors
import numpy as np
import astropy
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.optimize import curve_fit
from scipy import stats

from lmfit import (Minimizer, Parameters, conf_interval, conf_interval2d,
                   report_ci, report_fit)


# constant for transformation from kpc/Myr to km/s
trafo = 1000 / 1.023  # 1 kpc/Myr = (1000 / 1.023) km/s

# scaling factor for 3D plot
scale = 3

# age of the universe
t0 = 13800

# sun velocity vector in galactic rest frame
sun_vel = [11.1, 12.24 + 240, 7.25]

# Andromeda position in galactic coordinates
andromeda_pos = [121.172, -21.573]

andromeda_gal = SkyCoord(
    andromeda_pos[0], andromeda_pos[1], frame="galactic", unit="deg")


#
'''READ-OUT DATA FILES: '''
#

# readout process nr 1
# Read out the data of all distances at t_0
# readout1 = np.loadtxt('all-min-distances-to-sun.txt', dtype=str)
readout1 = np.loadtxt('/home/lguelzow/Nextcloud/MA/MA_Paper/HPC_Results/Plotting/present-time-half-mass.txt', dtype=str)

# Convert the whole table into an array of strings
dist1 = np.array(readout1)

# Get the travel time from the arrival time and send-off time
travel_time_present_time_half = [float(dist1[x, 2]) - float(dist1[x, 3]) for x in range(len(dist1))]


# readout process nr 2
# Read out the data of all distances at minimum distance
readout2 = np.loadtxt('/home/lguelzow/Nextcloud/MA/MA_Paper/HPC_Results/Plotting/minimum-distances-half-mass.txt', dtype=str)

# Convert the whole table into an array of strings
dist2 = np.array(readout2)

# Get the travel time from the arrival time and send-off time
travel_time_min_half = [float(dist2[x, 2]) - float(dist2[x, 3]) for x in range(len(dist2))]
# print(len(distance))
# print(len(MW_centre))


# readout process nr 3
# Read out the data of all distances at t_0
# readout1 = np.loadtxt('all-min-distances-to-sun.txt', dtype=str)
readout3 = np.loadtxt('/home/lguelzow/Nextcloud/MA/MA_Paper/HPC_Results/Plotting/present-time-same-mass.txt', dtype=str)

# Convert the whole table into an array of strings
dist3 = np.array(readout3)

# Get the travel time from the arrival time and send-off time
travel_time_present_time_same = [float(dist3[x, 2]) - float(dist3[x, 3]) for x in range(len(dist3))]


# readout process nr 4
# Read out the data of all distances at minimum distance
readout4 = np.loadtxt('/home/lguelzow/Nextcloud/MA/MA_Paper/HPC_Results/Plotting/minimum-distances-same-mass.txt', dtype=str)

# Convert the whole table into an array of strings
dist4 = np.array(readout4)

# Get the travel time from the arrival time and send-off time
travel_time_min_same = [float(dist4[x, 2]) - float(dist4[x, 3]) for x in range(len(dist4))]

# print(len(distance))
# print(len(MW_centre))

# Define confidence interval.
ci = 0.9973
# Convert to percentile point of the normal distribution.
# See: https://en.wikipedia.org/wiki/Standard_score
pp = (1. + ci) / 2.
# Convert to number of standard deviations.
nstd = 3 # stats.norm.ppf(pp)
print('Number of standard deviations', nstd)



#
# PLOTTING
# 

b = 20
b = np.arange(900, 3750, 150)

plt.figure(1)
# plot histogram of travel times
plt.hist(travel_time_present_time_same, bins=b, label='Equal-mass scenario', histtype="step", linewidth=2, color='black')
plt.hist(travel_time_present_time_half, bins=b, label='Half-mass scenario', histtype="step", linewidth=2, color='black', alpha=0.5)

# axis labels
plt.xlabel('Flight time from Andromeda to Milky Way [Myr]')
plt.ylabel('# of HVSs')

# title for plot
# plt.title("Flight times of present time HVSs")
plt.legend()
plt.savefig('Flight_times.pdf')
plt.show()
plt.close()

# plt.ylim(0, 77)
# plt.yscale('log')
# plt.xlim(0, 33)


plt.figure(2)
# plot histogram of travel times
plt.hist(travel_time_min_same, bins=b, label='Equal-mass scenario', histtype="step")
plt.hist(travel_time_min_half, bins=b, label='Half-mass scenario', histtype="step")

# axis labels
plt.xlabel('Travel time from Andromeda to Milky Way [Myr]')
plt.ylabel('# of HVSs')

# title for plot
plt.title("HVSs at minimum distance")
plt.legend()

# plt.ylim(0, 77)
# plt.yscale('log')
# plt.xlim(0, 33)

plt.show()