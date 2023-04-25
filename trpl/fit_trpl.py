"""
fit_trpl.py

Brief: Plot TRPL data for a given temperature and perform biexponential fit for lowest power.

Usage: python fit_trpl.py filename.csv (window) (gamma)

Author: Ian S. Dunn, Columbia University
"""

# Import modules.
import numpy as np
import sys
import matplotlib.pyplot as plt
import bottleneck as bn
import argparse
from util import *
from import_trpl import *
from model import *
import matplotlib as mpl
import pandas as pd
#import plotsettings

# Parse command line arguments.
parser = argparse.ArgumentParser(description='Process TRPL data.')
parser.add_argument('csv_file', type=str)
parser.add_argument('-w', '--window', type=int, default=-1)
parser.add_argument('-e', '--exclude', type=int, default=0)
parser.add_argument('-g', '--gamma', type=float, default=-1.)
parser.add_argument('-a', '--norm_average', action='store_true')
parser.add_argument('-l', '--log', action='store_true')
parser.add_argument('-p', '--plot_only', action='store_true')
parser.add_argument('-P', '--powers', action='store_true')
parser.add_argument('-D', '--densities', action='store_true')
parser.add_argument('-n', '--segments', default=1)
parser.add_argument('-m', '--maxtime', type=float, default=1000)
parser.add_argument('-r', '--margin', type=float, default=0.1)
parser.add_argument('-s', '--density_scaling_factor', type=float, default=1.54e14)
parser.add_argument('-x','--suppress_excimer_annihilation', action='store_true')
parser.add_argument('-f', '--fitting_power', type=int, default=0)

args = parser.parse_args()

# TRPL data filename.
csv_file = args.csv_file

# Record temperature.
T = int(args.csv_file.split('.')[0])

# Smoothing window.
window = args.window

# Single power to fit annihilation rate against.
fitting_power = args.fitting_power

# Number of powers to exclude from gamma fitting.
exclude = args.exclude
if fitting_power:
    exclude = 0

# Number of gamma fitting segments.
segments = int(args.segments)

log = args.log
plot_only = args.plot_only
norm_average = args.norm_average

if args.powers is True or args.densities is False:
    units = 'powers'
else:
    units = 'densities'

# Import normalized TRPL data and time in ns.
t, trpl = import_trpl(csv_file, norm_average)

dt = t[1] - t[0]
assert dt > 0

# Number of powers.
N = trpl.shape[0]

# Array to holds smoothed TRPL data.
smoothed = np.zeros_like(trpl)

# Load powers.
power_file = args.csv_file.split('.')[0]+'_' + units + '.txt'
powers = np.loadtxt(power_file)
#powers=np.delete(powers,1)
#print(powers)

# Scale powers from nW to cm^-3.
if units == 'powers':
    powers *= args.density_scaling_factor

# Set annihilation of slower component (excimer) to zero.
suppress_excimer_annihilation=args.suppress_excimer_annihilation
    
# Smooth data and plot.
for i in range(N):

    #smoothed[i] = bn.move_mean(trpl[i], window=window, min_count=1)
    from scipy.signal import savgol_filter

    if window > 1:
        smoothed[i] = savgol_filter(trpl[i], window_length=window, polyorder=1)
        smoothed[i] /= smoothed[i][0]
    else:
        smoothed[i] = bn.move_mean(trpl[i], window=1, min_count=1)
        #smoothed[i] = trpl[i]
        #smoothed[i] /= smoothed[i][0]

#std = bn.move_std(trpl - smoothed, window=window, min_count=1)
#std = bn.move_mean(std, window=50, min_count=1)
#std = std[:, [20]]

for i in range(N - exclude):

    where = np.where(t <= 20)

    if log is True:
        plt.semilogy(t[where], (smoothed[i][where]), color=plt.cm.plasma(float(i + 0.5) / N), label='$' + latex_float(powers[i]) + '\;cm^{-3}$')
    else:
        plt.plot(t, (smoothed[i]), color=plt.cm.plasma(float(i + 0.5) / N), label='$' + "{:.2e}".format(powers[i]) + '\;cm^{-3}$')

plt.ylabel('Normalized intensity', fontsize=16)
plt.xlabel('Time (ns)', fontsize=16)
plt.xlim(0, 10)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.text(0.6, 0.016, 'T='+str(T) +'K',fontsize=16)
plt.legend(fontsize=14)
#plt.savefig("%003d.pdf" % T, bbox='tight', dpi=300)
plt.show()

fig=plt.figure(figsize=(3.5,2.2))
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 7
plt.rcParams['axes.linewidth'] = 0.5

if plot_only:
    sys.exit()

# Fit low-power data to biexponential.
popt, pcov = curve_fit(biexponential, t, smoothed[0])

# Fitting parameters.
a, b, c, d = popt
b *= -1
d *= -1

# Order exponentials by lifetime.
if b < d:
    tmp = a
    a = c
    c = tmp
    tmp = b
    b = d
    d = tmp

# Print parameters.
print("k1 = %3.2f, error = %3.3f ns" % (b, np.sqrt(pcov[1,1])))
print("k2 = %3.2f, error = %3.3f ns" % (d, np.sqrt(pcov[3,3])))
print("A1 = %3.2f, t1 = %3.2f ns" % (a, 1. / b))
print("A2 = %3.2f, t2 = %3.2f ns" % (c, 1. / d))

# Plot biexponential fitting.
plt.plot(t, smoothed[0], color='b', label='TRPL data')
plt.plot(t, a * np.exp(-b * t), 'r', label='$\\tau_1='+'%.3f'%(1./b)+'$ ns')
plt.plot(t, c * np.exp(-d * t), 'k', label='$\\tau_2='+'%.3f'%(1./d)+'$ ns')
plt.plot(t, biexponential(t, a, -b, c, -d), color='g', label='Biexponential')
plt.title('Low-power biexponential fit to TRPL, T = ' + str(T) + ' K')
plt.xlabel('Time (ns)')
plt.ylabel('Normalized intensity')
plt.legend()
plt.tight_layout()
plt.show()

############ Section of code that deals with annihilation follows. ##############

try:
    assert powers.ndim == 1
    assert len(powers) == N
except:
    raise IOError(power_file + ' does not have the correct format.')

# Normalize all powers or densities to lowest power or density.
n0_relative = powers / powers[0]

truncation = np.where(t < args.maxtime)

t_trunc = t[truncation[0]]
smoothed_trunc = smoothed[:, truncation[0]]
#std = std[:, truncation[0]]


#fitting second highest power


def fit(segments):

    assert segments == 1 or len(t) == len(t_trunc)

    # Break time into segments.
    t_segments = time_segments(t_trunc, segments)

    # Starting exciton densities at each time segment.
    n_start = np.array([a, c])

    #if fitting_power == 0:
    n_start = np.tile(n_start, (segments + 1, N - exclude, 1)) * n0_relative[:N - exclude, np.newaxis]
    #else:
    #    n_start = np.tile(n_start, (segments + 1, 2, 1)) * n0_relative[[0, fitting_power], np.newaxis]

    # Input annihilation constant.
    gamma = args.gamma
     
    model = Model(b, d, n_start, suppress_excimer_annihilation)

    # Else calculate annihilation constant.
    if gamma < 0.:

        if fitting_power == 0:
            smoothed_relative = (smoothed_trunc * n0_relative[:, np.newaxis])[:N-exclude]
        else:
            smoothed_relative = (smoothed_trunc * n0_relative[:, np.newaxis])[[0, fitting_power]]

        smoothed_relative_segments = segments_like(t_segments, smoothed_relative)

        model.fitting_power = fitting_power
        gamma = model.fit_scaled_gamma(t_segments, smoothed_relative_segments) / powers[0]
        model.fitting_power = 0

        print('gamma = ', gamma[0] * 1e9)
        print('gamma error = ', gamma[1] * 1e9)

    else:
        gamma = np.array([gamma]) * 1e-9

    #std_trunc = (std * n0_relative[:, np.newaxis])[:N-exclude]
    #res = model.maximize_scaled_gamma(t_trunc, smoothed_relative, args.margin) / powers[0]
    #res = np.array([res])

    # Evaluate model.
    model_data = model.trpl_model(t_segments, gamma * powers[0]) / n0_relative[:, np.newaxis][:N-exclude]

    return gamma * 1e9, model_data, t_segments

gamma, model_data1, t_segments = fit(segments)

#print(model_data[3])
fig=plt.figure(figsize=(3.5,2.2))
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 6
plt.rcParams['axes.linewidth'] = 0.5

plt.plot(t, trpl[0], color=plt.cm.plasma(float(0 + 0.5) / N), label='$' + latex_float(powers[0]) + '\;cm^{-3}$',alpha=0.5)
plt.plot(t_trunc, model_data1[0], 'k:')
plt.plot(t, trpl[1], color=plt.cm.plasma(float(2 + 0.5) / N), label='$' + latex_float(powers[1]) + '\;cm^{-3}$',alpha=0.5)
plt.plot(t_trunc, model_data1[1], 'k:')


#fitting highest power

def fit2(segments):

    assert segments == 1 or len(t) == len(t_trunc)

    # Break time into segments.
    t_segments = time_segments(t_trunc, segments)

    # Starting exciton densities at each time segment.
    n_start = np.array([a, c])
    n_start = np.tile(n_start, (segments + 1, N - exclude, 1)) * n0_relative[:N - exclude, np.newaxis]

    # Input annihilation constant.
    gamma = args.gamma

    model = Model(b, d, n_start,suppress_excimer_annihilation)

    # Else calculate annihilation constant.
    if gamma < 0.:

        smoothed_relative = (smoothed_trunc * n0_relative[:, np.newaxis])[:N-exclude]

        smoothed_relative_segments = segments_like(t_segments, smoothed_relative)

        gamma = model.fit_scaled_gamma(t_segments, smoothed_relative_segments) / powers[0]

        print(gamma * 1e9)

    else:
        gamma = np.array([gamma]) * 1e-9

    #std_trunc = (std * n0_relative[:, np.newaxis])[:N-exclude]
    #res = model.maximize_scaled_gamma(t_trunc, smoothed_relative, args.margin) / powers[0]
    #res = np.array([res])

    # Evaluate model.
    model_data = model.trpl_model(t_segments, gamma * powers[0]) / n0_relative[:, np.newaxis][:N-exclude]

    return gamma * 1e9, model_data, t_segments

gamma, model_data2, t_segments = fit2(segments)


# Plot simulated power-dependent TRPL with measured power-dependent TRPL.

plt.plot(t, trpl[2], color=plt.cm.plasma(float(4 + 0.5) / N), label='$' + latex_float(powers[2]) + '\;cm^{-3}$',alpha=0.5)

plt.plot(t_trunc, model_data2[2], 'k:')

plt.ylabel('Normalized intensity')
plt.xlabel('Time (ns)')
plt.xlim(-0.1, 4)
plt.xticks(np.arange(0,4.1,1))
plt.yticks(np.arange(0,1.1,0.5))

#plt.ylim(10**(),1.1)
plt.legend(loc=1)
plt.tight_layout()
plt.savefig('fittedgamma_T'+str(T) +'.pdf', dpi=1000, bbox_inches='tight')
plt.show()




