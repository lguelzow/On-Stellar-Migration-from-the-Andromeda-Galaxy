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
from lmfit import (Minimizer, Parameters, conf_interval, conf_interval2d,
                   report_ci, report_fit)


# readout of main result file
# Read out the data from the result file
# readout = np.loadtxt('smallest_min_distances_total.txt', dtype=np.str)
# readout = np.loadtxt('Marchetti_v_data.txt', dtype=np.str)
readout = np.loadtxt('/home/lguelzow/Nextcloud/MA/MA_Paper/Paper_calcs/Marchetti-data/Marchetti_v_data.txt', dtype=str)

# Convert the whole table into an array of strings
result = np.array(readout)

# modify this later for sorted histogram by era and velocity
# Define the velocity magnitude array while excluding the 1st line
list_length = len(result) - 1


# Get the distance values from the read-out table and write into list
MW_vel = np.array([float(result[i + 1, 0]) for i in range(list_length)])

'''
# filter whole list for the stars where the exponential curve begins
MW_vel_fast = [MW_vel[i] for i in range(list_length) if MW_vel[i] >= 350]

# print(len(MW_vel_fast))


# load data for histogram fit
data = np.loadtxt("/home/lguelzow/Nextcloud/MA/MA_Paper/Paper_calcs/Marchetti-data/dataMarchetti_clipped.d", unpack=True)

MW_vel_fast = data[0]

# velocity cutoffs
v_cut_low = 400
'''

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
# readout1 = np.loadtxt('minimum-distances-same-mass.txt', dtype=str)
readout1 = np.loadtxt('/home/lguelzow/Nextcloud/MA/MA_Paper/HPC_Results/Plotting/present-time-same-mass.txt', dtype=str)

# Convert the whole table into an array of strings
dist1 = np.array(readout1)

# list of present velocity of HVS in same mass scenario
present_velocity_magnitude_same_mass_Gal_rest = np.array([trafo * np.sqrt((float(dist1[x, 9]) + sun_vel[0] / trafo) ** 2 \
            + (float(dist1[x, 10]) + sun_vel[1] / trafo) ** 2 + (float(dist1[x, 11]) + sun_vel[2] / trafo) ** 2) for x in range(len(dist1))])

present_velocity_magnitude_same_mass = [trafo * np.sqrt(float(dist1[x, 9]) ** 2 + float(dist1[x, 10]) ** 2 + float(dist1[x, 11]) ** 2) for x in range(len(dist1))]

# readout initial send-off times
sendoff_time_same_mass = np.array([float(dist1[x, 3]) for x in range(len(dist1))])

# readout process nr 2
# Read out the data of all distances at minimum distance
# readout2 = np.loadtxt('minimum-distances-half-mass.txt', dtype=str)
readout2 = np.loadtxt('/home/lguelzow/Nextcloud/MA/MA_Paper/HPC_Results/Plotting/present-time-half-mass.txt', dtype=str)

# Convert the whole table into an array of strings
dist2 = np.array(readout2)

# list of present velocity oh HVS in half mass scenario
present_velocity_magnitude_half_mass_Gal_rest = np.array([trafo * np.sqrt((float(dist2[x, 9]) + sun_vel[0] / trafo) ** 2 + (float(dist2[x, 10]) + sun_vel[1] / trafo) ** 2 + (float(dist2[x, 11]) + sun_vel[2] / trafo) ** 2) for x in range(len(dist2))])

present_velocity_magnitude_half_mass = np.array([trafo * np.sqrt(float(dist2[x, 9]) ** 2 + float(dist2[x, 10]) ** 2 + float(dist2[x, 11]) ** 2) for x in range(len(dist2))])

# read out intital sendoff times
sendoff_time_half_mass = np.array([float(dist2[x, 3]) for x in range(len(dist2))])

# print(np.amin(sendoff_time_half_mass))

# print(len(present_velocity_magnitude_same_mass))
# print(len(sendoff_time_half_mass))
# print(len(sendoff_time_same_mass))

# velocity cutoffs
mask = (present_velocity_magnitude_half_mass_Gal_rest > 400) & (present_velocity_magnitude_half_mass_Gal_rest < 900)
mask2 = present_velocity_magnitude_same_mass_Gal_rest > 400 & (present_velocity_magnitude_same_mass_Gal_rest < 900)
v_cut_low = 400

mask_sendoff_same = (sendoff_time_same_mass < 11000) & (sendoff_time_same_mass > 10000)
mask_sendoff_half = (sendoff_time_half_mass < 11000) & (sendoff_time_half_mass > 10000)

# print(len(present_velocity_magnitude_same_mass_Gal_rest[mask_sendoff_same]))
# print(len(present_velocity_magnitude_half_mass_Gal_rest[mask_sendoff_half]))


plt.figure(1)
# plt.hist(test_plot, 30)
# plt.hist(present_velocity_magnitude_same_mass_Gal_rest[mask2], bins=5, histtype='step', linewidth=2, color='same', label='Equal-mass scenario')
plt.xlabel('Velocity in Galactic rest frame [km/s]')
plt.ylabel('# of HVSs')
plt.title('HVSs at present time')
# plt.xlim(-50, 1000)
# plt.ylim(0.3, 5 * 10 ** 6)
# plt.figure(2)
# plt.hist(test_plot, 30)

# fit stuff

# plt.figure(2)
# plt.hist(test_plot, 30)

# bin number
b = 8
b_same = np.arange(375, 1125, 75)
b_half = np.arange(375, 1125, 75)

b_same_fit = np.arange(450, 975, 75)
b_half_fit = np.arange(450, 975, 75)

bin_marchetti = np.arange(300, 1100, 75)


# 400 low cut version
# b_same = np.arange(375, 1125, 75)
# b_half = np.arange(375, 1125, 75)

# b_same_fit = np.arange(375, 975, 75)
# b_half_fit = np.arange(375, 975, 75)

# bin_marchetti = np.arange(300, 1100, 75)

# indices determining confidence levels
# 1st index determines confidence level (0: 99.73%; 1: 95.45%; 2: 68.27%, 3: BEST VALUE)
# 1 - 3 are the low bounds and 4 - 6 are the high bounds (reverse order)
index_high = 5
index_low = 1

# give bins as a linspace
# b = np.arange(400, 1100, 75)

# same fit for other data set
n, bins, _ = plt.hist(present_velocity_magnitude_same_mass_Gal_rest[mask2], \
                       bins=b_same_fit, histtype='stepfilled', linewidth=2, color='blue', \
                        label='Equal-mass scenario', alpha=0.35)

# correction so that x and y for the fit have the same shape
# take value in the middle of bins instead of edges
xdata = [bins[i] + ((bins[1] - bins[0]) / 2) for i in range(len(bins))]
del xdata[-1]
xdata_same = np.array(xdata)
ydata = np.array(n)

print(xdata_same)
print(n)

# least square fit; fit function
def residual(params, v, data):
    sigma_v = params['sigma_v']
    A = params['amplitude']
    model = (A / sigma_v) * np.exp((-1) * ((v - v_cut_low) / sigma_v))
    return (data - model) / np.sqrt(data)

# calculate the chisq for the model p = 1/d
params = Parameters()
params.add('sigma_v', value=59.7846854, vary=True)
params.add('amplitude', value=94778.1540, vary=True)

mini = Minimizer(residual, params, fcn_args=(xdata_same, ydata))
result0 = mini.leastsq()

fit_res0 = residual(result0.params, xdata_same, ydata)
report_fit(result0)

ci = conf_interval(mini, result0)
print(ci)
report_ci(ci)

# actual fit data that can be plotted
fit_same = ((result0.params['amplitude'].value / result0.params['sigma_v'].value) * np.exp((-1) \
                         * ((xdata_same - v_cut_low) / result0.params['sigma_v'].value)))
# plot the fit
plt.plot(xdata_same, fit_same, '--', color='red',
         label='Exponential fit', linewidth=2.5)  # , alpha=0.5)

# define confidence interval functions
# get conidence interval values from ci:
# 1st index determines confidence level (0: 99.73%; 1: 95.45%; 2: 68.27%, 3: BEST VALUE)
# 1 - 3 are the low bounds and 4 - 6 are the high bounds (reverse order)
# 2nd index gives you the value of the confidence level

ci_high_same = ci["amplitude"][index_high][1] / ci["sigma_v"][index_high][1] \
           * np.exp((-1) * ((xdata_same - v_cut_low) / ci["sigma_v"][index_high][1]))
ci_low_same = ci["amplitude"][index_low][1] / ci["sigma_v"][index_low][1] \
           * np.exp((-1) * ((xdata_same - v_cut_low) / ci["sigma_v"][index_low][1]))

# plot confidence regions
plt.fill_between(xdata_same, ci_low_same, ci_high_same, interpolate=True, color='blue', alpha=0.1)



# same fit for the other data set
n, bins, _ = plt.hist(present_velocity_magnitude_half_mass_Gal_rest[mask], \
                      bins=b_half_fit, histtype='stepfilled', linewidth=2, color='black', \
                        label='Half-mass scenario', alpha=0.35)

# correction so that x and y for the fit have the same shape
# take value in the middle of bins instead of edges
xdata = [bins[i] + ((bins[1] - bins[0]) / 2) for i in range(len(bins))]
del xdata[-1]
xdata_half = np.array(xdata)
ydata = np.array(n)

print(xdata_half)
print(ydata)
print(n)

# least square fit; fit function
def residual(params, v, data):
    sigma_v = params['sigma_v']
    A = params['amplitude']
    model = (A / sigma_v) * np.exp((-1) * ((v - v_cut_low) / sigma_v))
    return (data - model) / np.sqrt(data)

# calculate the chisq for the model p = 1/d
params = Parameters()
params.add('sigma_v', value=59.7846854, vary=True)
params.add('amplitude', value=94778.1540, vary=True)

mini = Minimizer(residual, params, fcn_args=(xdata_half, ydata))
result0 = mini.leastsq()

fit_res0 = residual(result0.params, xdata_half, ydata)
report_fit(result0)

ci = conf_interval(mini, result0)
report_ci(ci)


# actual fit data that can be plotted
fit_half = ((result0.params['amplitude'].value / result0.params['sigma_v'].value) * np.exp((-1) \
             * ((xdata_half - v_cut_low) / result0.params['sigma_v'].value)))
# plot the fit
plt.plot(xdata_half, fit_half, '--', color='black',
         label='Exponential fit', linewidth=2.5)  # , alpha=0.5)

# define confidence interval functions
# get conidence interval values from ci:
# 1st index determines confidence level (0: 99.73%; 1: 95.45%; 2: 68.27%, 3: BEST VALUE)
# 1 - 3 are the low bounds and 4 - 6 are the high bounds (reverse order)
# 2nd index gives you the value of the confidence level instead of the fit
ci_high_half = ci["amplitude"][index_high][1] / ci["sigma_v"][index_high][1] \
           * np.exp((-1) * ((xdata_half - v_cut_low) / ci["sigma_v"][index_high][1]))
ci_low_half = ci["amplitude"][index_low][1] / ci["sigma_v"][index_low][1] \
           * np.exp((-1) * ((xdata_half - v_cut_low) / ci["sigma_v"][index_low][1]))

# plot confidence regions
plt.fill_between(xdata_half, ci_low_half, ci_high_half, interpolate=True, color='black', alpha=0.1)

plt.yscale('log')
plt.legend()
plt.savefig('Velocity-dist-Gal-restframe-fit.pdf')
plt.show()
plt.close()




plt.figure(2)
# plt.hist(test_plot, 30)
plt.hist(present_velocity_magnitude_same_mass_Gal_rest, bins=b_same, histtype='step', \
          linewidth=1.5, color='black', label='Equal-mass scenario', alpha=1)

# plt.hist(present_velocity_magnitude_same_mass_Gal_rest, bins=35, histtype='stepfilled', \
#           linewidth=1.5, color='blue', label='test_same', alpha=0.5)

plt.xlabel(r'Velocity in Galactic rest frame [kms$^{-1}$]')
plt.ylabel('# of HVSs')
plt.yscale('log')
# plt.xlim(360, 715)
plt.xlim(340, 1080)
# plt.ylim(1.5, 800)
plt.ylim(0.2, 100000)

plt.hist(present_velocity_magnitude_half_mass_Gal_rest, bins=b_half, histtype='step', \
          linewidth=1.5, color='black', label='Half-mass scenario', alpha=0.5)

# test histograms to show the data in finer bins
'''
plt.hist(present_velocity_magnitude_half_mass_Gal_rest, bins=35, histtype='stepfilled', \
          linewidth=1.5, color='black', label='test_half', alpha=0.5)
plt.hist(present_velocity_magnitude_same_mass_Gal_rest, bins=35, histtype='stepfilled', \
          linewidth=1.5, color='blue', label='test_same', alpha=0.5)
'''

# insert fit functions
plt.plot(xdata_same, fit_same, '--', color='blue', linewidth=2)  # , alpha=0.5)

plt.plot(xdata_half, fit_half, '--', color='blue', label='Exponential fits', linewidth=2)  # , alpha=0.5)

# plot confidence regions
plt.fill_between(xdata_same, ci_low_same, ci_high_same, interpolate=True, color='blue', alpha=0.15, label=r'$2\sigma$ confidence regions')
# plot confidence regions
plt.fill_between(xdata_half, ci_low_half, ci_high_half, interpolate=True, color='blue', alpha=0.15)


# also add the Marchetti plot
mask_marchetti = (MW_vel > 300) & (MW_vel < 1050)

# find data by generating the histogram
plt.hist(MW_vel[mask_marchetti], bins=bin_marchetti, label='MW velocity distribution (see Fig. 2)',
                      color='black', histtype='step', linewidth=1, linestyle=':')

# plt.title("Only HVSs sent between 11 and 11.5 Gyrs")
plt.legend(loc='upper right')
plt.savefig('Velocity-dist-Gal-restframe-complete.pdf')


plt.figure(3)
# plt.hist(test_plot, 30)
plt.hist(present_velocity_magnitude_same_mass_Gal_rest[mask_sendoff_same], bins=np.arange(400, 800, 6.25), histtype='step', \
          linewidth=2, color='blue', label='Equal-mass scenario')
plt.xlabel('Velocity in LSR [km/s]')
plt.ylabel('# of HVSs')
# plt.yscale('log')
# plt.xlim(-50, 1000)
# plt.ylim(0.3, 5 * 10 ** 6)

plt.hist(present_velocity_magnitude_half_mass_Gal_rest[mask_sendoff_half], bins=np.arange(400, 800, 6.25), histtype='stepfilled', \
          linewidth=2, color='black', label='Half-mass scenario', alpha=0.5)

plt.title("Only HVSs sent off between 11 and 12 Gyrs")
plt.legend()
plt.savefig('Velocity-dist-Gal-restframe.pdf')
plt.show()



