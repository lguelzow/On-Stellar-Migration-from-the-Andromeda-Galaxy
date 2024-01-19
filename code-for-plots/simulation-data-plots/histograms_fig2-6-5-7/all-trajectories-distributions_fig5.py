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

# readout process nr 1
# Read out the data of all distances at t_0
# readout1 = np.loadtxt('/home/lguelzow/Nextcloud/MA/MA_Paper/Paper_github/On-Stellar-Migration-from-the-Andromeda-Galaxy/simulation-results/power_law_test/combined/all-minimum-distances-power-law.txt', dtype=np.str)
readout1 = np.loadtxt('/home/lguelzow/Nextcloud/MA/MA_Paper/Paper_github/On-Stellar-Migration-from-the-Andromeda-Galaxy/simulation-results/no_plummer_test/all-minimum-distances-no-plummer.txt', dtype=np.str)

# readout1 = np.loadtxt('all-min-distances-same-mass.txt', dtype=np.str)
#  readout2 = np.loadtxt('all-min-distances-half-mass.txt', dtype=np.str)

# Convert the whole table into an array of strings
dist1 = np.array(readout1)
# dist2 = np.array(readout2)

# Define the bigger velocity array while excluding the 1st line
distance_same = np.array([float(dist1[x]) for x in range(len(dist1))])
# distance_half = np.array([float(dist2[x]) for x in range(len(dist2))])

# list for stars within 50kpc of Milky Way centre
MW_centre_same = [distance_same[i] for i in range(len(dist1)) if distance_same[i] < 50]

# MW_centre_half = [distance_half[i] for i in range(len(dist2)) if distance_half[i] < 50]

print("Same-mass scenario:")
print(len(distance_same))
print(len(MW_centre_same))

print("\n" + "Half-mass scenario")
# print(len(distance_half))
# print(len(MW_centre_half))


plt.figure(1)
# plt.hist(test_plot, 30)
plt.hist(distance_same, bins=np.arange(0, 1200, 40), histtype='step', color='black', linewidth=2, label="Equal-mass scenario")
# plt.hist(distance_half, bins=np.arange(0, 1200, 40), histtype='step', color='black', linewidth=2, label="Half-mass scenario", alpha=0.5)
plt.xlabel('Distance $r$ to the MW centre [kpc]')
plt.ylabel('# of HVSs')
plt.yscale('log')
plt.legend(loc="upper left")
plt.xlim(-30, 1180)
plt.savefig('SIM-All-stars-distance_histogram.pdf')
plt.show()

'''
plt.figure(2)
# plt.hist(test_plot, 30)
plt.hist(MW_centre, bins=25, histtype='step')
plt.xlabel('Distance $r$ to the MW centre [kpc]')
plt.ylabel('#of HVSs')
# plt.yscale('log')
# plt.xlim(0, 33)
plt.show()
'''
