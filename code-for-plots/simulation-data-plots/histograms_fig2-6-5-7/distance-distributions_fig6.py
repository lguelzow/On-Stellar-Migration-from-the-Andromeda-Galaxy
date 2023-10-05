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
# readout1 = np.loadtxt('all-min-distances-to-sun.txt', dtype=np.str)
readout1 = np.loadtxt('/home/lguelzow/Nextcloud/MA/MA_Paper/HPC_Results/Plotting/present-time-half-mass.txt', dtype=np.str)

# Convert the whole table into an array of strings
dist1 = np.array(readout1)

# Get the distance values from the read-out table and convert the strings into floats
present_distance_half = [float(dist1[x, 0]) for x in range(len(dist1))]


# readout process nr 2
# Read out the data of all distances at minimum distance
readout2 = np.loadtxt('/home/lguelzow/Nextcloud/MA/MA_Paper/HPC_Results/Plotting/minimum-distances-half-mass.txt', dtype=np.str)

# Convert the whole table into an array of strings
dist2 = np.array(readout2)

# Get the distance values from the read-out table and convert the strings into floats
min_distance_half = [float(dist2[x, 0]) for x in range(len(dist2))]
# print(len(distance))
# print(len(MW_centre))


# readout process nr 3
# Read out the data of all distances at t_0
# readout1 = np.loadtxt('all-min-distances-to-sun.txt', dtype=np.str)
readout3 = np.loadtxt('/home/lguelzow/Nextcloud/MA/MA_Paper/HPC_Results/Plotting/present-time-same-mass.txt', dtype=np.str)

# Convert the whole table into an array of strings
dist3 = np.array(readout3)

# Get the distance values from the read-out table and convert the strings into floats
present_distance_same = [float(dist3[x, 0]) for x in range(len(dist3))]


# readout process nr 4
# Read out the data of all distances at minimum distance
readout4 = np.loadtxt('/home/lguelzow/Nextcloud/MA/MA_Paper/HPC_Results/Plotting/minimum-distances-same-mass.txt', dtype=np.str)

# Convert the whole table into an array of strings
dist4 = np.array(readout4)

# Get the distance values from the read-out table and convert the strings into floats
min_distance_same = [float(dist4[x, 0]) for x in range(len(dist4))]

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

plt.figure(0)
b = 20
# plt.hist(test_plot, 30)
# ,  weights=np.ones(len(min_D)) / len(min_D))
n, bins, _ = plt.hist(present_distance_same, bins=b, label=None)


# define cubic function
def fit_func(x, a, b, c):
    # Curve fitting function
    return a * x ** 3 + b * x ** 2 + c * x  # d=0 is implied

# define cubic function
def fit_func_quad(x, b):
    # Curve fitting function
    return b * x ** 2#  + c * x  # d=0 is implied


# Curve fitting
params = curve_fit(fit_func, bins[0:b], n[0:b])
[a, b, c] = params[0]
x_fit = np.linspace(bins[0], bins[-1], 100)
y_fit = a * x_fit ** 3 + b * x_fit ** 2 + c * x_fit
# plot fit
plt.plot(x_fit, y_fit, c='orange', label='$r^3$ function fit')
plt.xlabel('Distance $r$ to the Sun [kpc]')
plt.ylabel('# of HVSs')
# plt.ylim(0, 77)
# plt.yscale('log')
# plt.xlim(0, 33)
# plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.show()


plt.figure(3)
# plt.hist(test_plot, 30)

'''n1, bins1, _ = plt.hist(present_distance_same, bins=np.arange(0, 50, 3+1/3), histtype='step', linewidth=2, color='black', label='Equal-mass scenario')
n2, bins2, _ = plt.hist(present_distance_half, bins=np.arange(0, 50, 3+1/3), histtype='step', linewidth=2, color='black', label='Half-mass scenario', alpha=0.5)

bins1 = np.array([bins1[i] + (bins1[1] - bins1[0]) / 2 for i in range(len(bins1) - 1)])
bins2 = np.array([bins2[i] + (bins2[1] - bins2[0]) / 2 for i in range(len(bins2) - 1)])

# define residual for fit evaluation
# least square fit; fit function
def residual(params, x, data):
    b = params['b_quad']
    model = b * x ** 2
    return (data - model) / np.sqrt(data)

params = Parameters()
params.add('b_quad', value=0.1, vary=True)

# calculate fits
mini1 = Minimizer(residual, params, fcn_args=(bins1, n1))
result1 = mini1.leastsq()

mini2 = Minimizer(residual, params, fcn_args=(bins2, n2))
result2 = mini2.leastsq()

# print fit results
fit_res1 = residual(result1.params, bins1, n1)
fit_res2 = residual(result2.params, bins2, n2)
report_fit(result1)
report_fit(result2)'''

'''# calculate confidence intervals
ci1 = conf_interval(mini1, result1)
print(ci1)
report_ci(ci1)

ci2 = conf_interval(mini2, result2)
print(ci2)
report_ci(ci2)'''


'''
# Curve fitting
params1 = curve_fit(fit_func, bins1[0:15], n1[0:15])
# print(params1)
params2 = curve_fit(fit_func, bins2[0:15], n2[0:15])
[a1, b1, c1] = params1[0]
[a2, b2, c2] = params2[0]
x_fit1 = np.linspace(bins1[0], bins1[-1], 100)
y_fit1 = a1 * x_fit1 ** 3 + b1 * x_fit1 ** 2 + c1 * x_fit1

x_fit2 = np.linspace(bins2[0], bins2[-1], 100)
y_fit2 = a2 * x_fit2 ** 3 + b2 * x_fit2 ** 2 + c2 * x_fit2
# plot fit
plt.plot(x_fit1, y_fit1, '--', c='grey',  label='Polynomial fits', linewidth=3)
plt.plot(x_fit2, y_fit2, '--',  c='grey', linewidth=3)
'''

plt.xlabel('Distance $r$ to the MW centre at $t_0$ [kpc]')
plt.ylabel('# of HVSs')
# plt.yscale('log')
plt.xlim(0, 52)
plt.legend(loc='upper left')
# plt.savefig('Present-distance-histogram.pdf')


# plt.histtest_plot, 30)
n1, bins1, _ = plt.hist(present_distance_same, bins=np.arange(0, 50, 3+1/3), histtype='step', linewidth=2, color='black', label='Equal-mass scenario')
n2, bins2, _ = plt.hist(present_distance_half, bins=np.arange(0, 50, 3+1/3), histtype='step', linewidth=2, color='black', label='Half-mass scenario', alpha=0.5)

bins1 = [bins1[i] + (bins1[1] - bins1[0]) / 2 for i in range(len(bins1) - 1)]
bins2 = [bins2[i] + (bins2[1] - bins2[0]) / 2 for i in range(len(bins2) - 1)]

# Curve fitting
params1 = curve_fit(fit_func_quad, bins1[0:15], n1[0:15], full_output=True)
print(params1)
params2 = curve_fit(fit_func_quad, bins2[0:15], n2[0:15], full_output=True)
print(params2)
b1 = params1[0] 
# b1 = result1.params['b_quad'].value
b2 = params2[0] 
# b2 = result2.params['b_quad'].value
x_fit1 = np.linspace(bins1[0], bins1[-1], 100)
y_fit1 = b1 * x_fit1 ** 2

x_fit2 = np.linspace(bins2[0], bins2[-1], 100)
y_fit2 = b2 * x_fit2 ** 2

# Standard deviation errors on the parameters.
perr1 = np.sqrt(np.diag(params1[1]))
perr2 = np.sqrt(np.diag(params2[1]))
# Add/subtract nstd standard deviations to parameters to obtain the upper/lower confidence interval.
popt_up1 = params1[0] + nstd * perr1
popt_dw1 = params1[0] - nstd * perr1

popt_up2 = params2[0] + nstd * perr2
popt_dw2 = params2[0] - nstd * perr2

y_fit1_up = popt_up1[0] * x_fit1 ** 2
y_fit1_down = popt_dw1[0] * x_fit1 ** 2

y_fit2_up = popt_up2[0] * x_fit2 ** 2
y_fit2_down = popt_dw2[0] * x_fit2 ** 2

# plot fit
plt.plot(x_fit1, y_fit1, '--', c='blue',  label='Quadratic fits', linewidth=2)
plt.plot(x_fit2, y_fit2, '--',  c='blue', linewidth=2)
plt.fill_between(x_fit1, y_fit1_up, y_fit1_down, interpolate=False, color='blue', alpha=0.15, label=r'$3\sigma$ confidence regions')
plt.fill_between(x_fit2, y_fit2_up, y_fit2_down, interpolate=False, color='blue', alpha=0.15)

plt.xlabel('Distance $r$ to the MW centre at $t_0$ [kpc]')
plt.ylabel('# of HVSs')
# plt.yscale('log')
# plt.ylim(0, 100)
plt.legend(loc='upper left')
plt.savefig('Present-distance-histogram-fit.pdf')


plt.figure(1)
# plt.hist(test_plot, 30)
plt.hist(present_distance_same, bins=np.arange(0, 50, 3+1/3), histtype='step',
         linewidth=2, color='black', label='Equal-mass scenario')
plt.hist(present_distance_half, bins=np.arange(0, 50, 3+1/3), histtype='step',
         linewidth=2, color='black', label='Half-mass scenario', alpha=0.5)
plt.xlabel('Distance $r$ to the MW centre at $t_0$ [kpc]')
plt.ylabel('# of HVSs')
# plt.yscale('log')
# plt.ylim(0, 100)
plt.legend(loc='upper left')
plt.savefig('Present-distance-histogram.pdf')


plt.figure(2)
# plt.hist(test_plot, 30)
plt.hist(min_distance_same, bins=18,
         histtype='step', linewidth=2, color='black', label='Equal-mass scenario')
plt.hist(min_distance_half, bins=18,
         histtype='step', linewidth=2, color='black', label='Half-mass scenario', alpha=0.5)
plt.xlabel('Minimum distance $r$ to the MW centre [kpc]')
plt.ylabel('# of HVSs')
# plt.yscale('log')
# plt.ylim(0, 110)
plt.legend(loc='upper left')
plt.savefig('Minimum-distance-histogram.pdf')
plt.show()
