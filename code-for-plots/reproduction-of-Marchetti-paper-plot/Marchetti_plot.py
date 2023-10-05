import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import astropy
from astropy.coordinates import SkyCoord
import astropy.units as u
import scipy.optimize
from lmfit import (Minimizer, Parameters, conf_interval, conf_interval2d,
                   report_ci, report_fit)


# readout of main result file
# Read out the data from the result file
# readout = np.loadtxt('smallest_min_distances_total.txt', dtype=np.str)
# readout = np.loadtxt('Marchetti_v_data.txt', dtype=np.str)
readout = np.loadtxt('Marchetti_v_data.txt', dtype=str)

# Convert the whole table into an array of strings
result = np.array(readout)

# modify this later for sorted histogram by era and velocity
# Define the velocity magnitude array while excluding the 1st line
list_length = len(result) - 1


# Get the distance values from the read-out table and write into list
MW_vel = [float(result[i + 1, 0]) for i in range(list_length)]


# filter whole list for the stars where the exponential curve begins
MW_vel_fast = [MW_vel[i] for i in range(list_length) if MW_vel[i] >= 350]

# print(len(MW_vel_fast))


# load data for histogram fit
data = np.loadtxt("dataMarchetti_clipped.d", unpack=True)

MW_vel_fast = data[0]

# velocity cutoffs
v_cut_low = 400


#
'''PLOTTING: '''
#

# HISTOGRAMS

# fast tail of velocity distribution
plt.figure(3)
# generate histogram for the data to fit
# number of bins
b = 20

# find data by generating the histogram
n, bins, _ = plt.hist(MW_vel_fast, bins=b, label='Star velocity distribution',
                      color='black', histtype='step', linewidth=1.5)

# print(n)
# print(len(bins))

xdata = [bins[i] + ((bins[1] - bins[0]) / 2) for i in range(len(bins))]

del xdata[-1]

xdata = np.array(xdata)
ydata = np.array(n)
print(xdata)
print(ydata)


# least square fit
def residual(params, v, data):
    sigma_v = params['sigma_v']
    A = params['amplitude']
    model = (A / sigma_v) * np.exp((-1) * ((v - v_cut_low) / sigma_v))
    return (data - model) / np.sqrt(data)


# calculate the chisq for the model p = 1/d
params = Parameters()
params.add('sigma_v', value=5.86910843e+01, vary=True)
params.add('amplitude', value=8.76217321e+04, vary=True)

mini = Minimizer(residual, params, fcn_args=(xdata, ydata))
result0 = mini.leastsq()

fit_res0 = residual(result0.params, xdata, ydata)
report_fit(result0)

ci = conf_interval(mini, result0)
# print(ci)
# print(ci["amplitude"][0][1])
report_ci(ci)

fit0 = ((result0.params['amplitude'].value / result0.params['sigma_v'].value) * np.exp((-1)
                                                                                * ((xdata - v_cut_low) / result0.params['sigma_v'].value)))
plt.plot(xdata, fit0, '--', color='red', linewidth=2, label='Exponential fit')

# define confidence interval functions
# get conidence interval values from ci:
# 1st index determines confidence level (0: 99.73%; 1: 95.45%; 2: 68.27%, 3: BEST VALUE)
# 1 - 3 are the low bounds and 4 - 6 are the high bounds (reverse order)
# 2nd index gives you the value of the confidence level
ci_high = ci["amplitude"][-1][1] / ci["sigma_v"][-1][1] \
           * np.exp((-1) * ((xdata - v_cut_low) / ci["sigma_v"][-1][1]))
ci_low = ci["amplitude"][0][1] / ci["sigma_v"][0][1] \
           * np.exp((-1) * ((xdata - v_cut_low) / ci["sigma_v"][0][1]))

# plot confidence regions
plt.fill_between(xdata, ci_low, ci_high)


# plot data and fits with a legend
plt.legend(loc='best')


plt.xlabel('Star velocity [km/s]')
plt.ylabel('# of stars')
# plt.yscale('log')
plt.xticks(np.arange(0, 900, step=100))
plt.xlim(360, 715)
plt.ylim(1.5, 800)
# plt.legend(loc='upper right', fontsize=8)
plt.savefig('Marchetti_fast_data2.pdf')


# distribution of star velocities from Marchetti data
plt.figure(2)
plt.rcParams['font.size'] = 10
# plt.hist(Marchetti_v_data_plot, 30)
plt.hist(MW_vel, 40, histtype='step', color='black',
         linewidth=1.5, label='Distribution of star velocities\n in the Milky Way', alpha=1)
# plot fit from previous plot into this histogram too
plt.plot(xdata, fit0, '--', color='blue',
         label='Exponential fit', linewidth=2)  # , sigma_v=0.5)
# plot confidence regions
plt.fill_between(xdata, ci_low, ci_high, color='blue', alpha=0.15, label=r'3$\sigma$ confidence region')
plt.xlabel(r'Velocity in Galactic rest frame [kms$^{-1}$]')
plt.ylabel('# of stars')
plt.yscale('log')
plt.xticks(np.arange(0, 900, step=200))
plt.xlim(-20, 900)
plt.legend(loc='upper right')
plt.savefig('Marchetti_v_data.pdf')
plt.show()

'''
# fast tail of velocity distribution(log)
plt.figure(2)
# plt.hist(Marchetti_v_data_plot, 30)
plt.hist(MW_vel_fast, 20)
plt.xlabel('Velocity [km/s]')
plt.ylabel('# of stars')
plt.yscale('log')
plt.xticks(np.arange(0, 900, step=100))
plt.xlim(350, 900)
plt.savefig('Marchetti_fast_data.pdf')


# fast tail of velocity distribution
plt.figure(3)
# plt.hist(Marchetti_v_data_plot, 30)
b = 28
# print(b)
n, bins, _ = plt.hist(MW_vel_fast, bins=b, label=None)

# bins = np.linspace(600.05134754, 880.78443527, 50)

# print(len(MW_vel_fast))
# print(n)

xs = bins

for i in range(len(bins)):
    xs[i] = bins[i] + ((bins[1] - bins[0]) / 2)

xs = xs.tolist()
del xs[-1]

xs = np.array(xs)
ys = np.array(n)

# print(xs)

# define exp function


def exponential(x, v0, y0, k):
    return y0 - (v0 / k) * (1 - np.exp(-k * x))


popt_exponential, pcov_exponential = scipy.optimize.curve_fit(
    exponential, xs, ys, p0=[10400, 633373, 0.016], maxfev=1000000)

# print(*popt_exponential)
# print(pcov_exponential)

# print(exponential(xs, *popt_exponential))
# print(exponential(900, 6000, 400000, 0.015))

# plt.plot(xs, ys, '.', label="Original data")
# plot exponential fit
plt.plot(xs, exponential(xs, *popt_exponential), label="Exponential fit")


plt.xlabel('Velocity [km/s]')
plt.ylabel('# of stars')
# plt.xscale('log')
plt.xticks(np.arange(0, 900, step=200))
plt.xlim(480, 900)
plt.ylim(0, 450)
plt.legend(loc='upper right', fontsize=8)
plt.savefig('Marchetti_fast_data2.pdf')


# test plot for fit function
plt.figure(4)

xtest = np.linspace(800, 1200, 1000)

# print(xtest)

plt.plot(xtest, exponential(xtest, *popt_exponential), label="Exponential fit")
plt.xticks(np.arange(800, 1200, step=50))
plt.xlim(800, 1150)
plt.ylim(0, 2)
plt.xlabel('Velocity [km/s]')
plt.ylabel('Amount of stars')
plt.legend(loc='upper right', fontsize=8)
plt.savefig('testplot.pdf')
'''
