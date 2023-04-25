import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import sys
import numpy as np
sys.path.append('..')
import Header
Colors, FSize, Vars = Header.InitFigure()
import pickle

cmap = matplotlib.cm.get_cmap('coolwarm_r')

from exciton import Exciton

fig = plt.figure(figsize=(Vars['DWidth'], 3))

gs = gridspec.GridSpec(1, 3, width_ratios=[2., 3., 1.], figure=fig)

def plot(molecule, periodic, color, index):

    ax = plt.subplot(gs[index])

    exciton = Exciton(molecule, 0.00595, periodic)

    #exciton_pickle = open(molecule + "_exciton.pickle", "wb")
    #import pickle
    #pickle.dump(exciton, exciton_pickle)
    
    #with open(molecule + "_exciton.pickle", "rb") as exciton_pickle:
    #    exciton = pickle.load(exciton_pickle)

    print('One exciton Hamiltonian eigen-energies:')
    print(exciton.w)

    vv = exciton.loaded_wavefunction()

    np.savetxt('projected_eigenvectors_' + molecule + '.txt', vv[:exciton.nsites * exciton.nvib][::exciton.nvib])

    sums = []

    colors = ['b'] * vv.shape[1]

    if molecule == 'n-ph-pdi':
        num_col = 2
    else:
        num_col = 3

    for i in range(vv.shape[1]):
        v = vv[:, i]
        #if exciton.w[i] * 1.4 > 300:
        #    continue
        sums.append(np.sum(v[:exciton.nsites * exciton.nvib][::exciton.nvib])**2)
        if sums[i] < 1.e-4:
            colors[i] = 'r'

    #sums = np.array(sums)
    sums /= np.max(sums)
    #sums=(sums-np.min(sums))/( np.max(sums)-np.min(sums))
    #sums = sums * 1.1
    #sums = np.log(sums) / np.log(1.1)

    for i in range(len(sums)):
        v = vv[:, i]
        ax.plot(i % num_col + 0.02 + np.arange(10.) / 10. * 0.83, v[:exciton.nsites * exciton.nvib][::exciton.nvib] * 20/1.4 + exciton.w[i], color=cmap(sums[i]), linewidth=2)
        #if sums[-1] > 1.e-4:
        #    if molecule == 'n-ph-pdi':
        #        ax.text(i - 0.05, exciton.w[i] * 1.4 + 23, "%0.3f" % sums[-1])
        #    else:
        #        ax.text(i - 0.05, exciton.w[i] * 1.4 + 23, "%0.3f" % sums[-1])

    ax.hlines(exciton.w, np.arange(len(exciton.w)) % num_col, np.arange(len(exciton.w)) % num_col + 0.8, 'k')

    if molecule == 'n-ph-pdi':
        plt.title('N-phenyl PDI')
    elif molecule == 'tp-pdi':
        plt.title('Tetraphenyl PDI')

    plt.ylim(-10/1.4, 300/1.4)
    plt.xticks([])

    if not index:
        plt.ylabel('$\\omega (cm^{-1})$')

    if index:
        plt.yticks([])

    return exciton.w, np.array(sums)

def rate(w, sums):

    T = np.arange(0.1, 300, 0.1)
    beta = 1./T
    boltzmann = np.exp(-beta * w[:, np.newaxis])
    boltzmann /= np.sum(boltzmann, axis=0)

    print("Relative exciton annihilation strength, relative exciton Boltzmann weight:")
    for i in range(len(sums)):
        print(sums[i], ",", boltzmann[i, -1])

    return T, np.dot(sums, boltzmann)

periodic = int(sys.argv[1])

w, s = plot('n-ph-pdi', periodic, 'b', 0)
T, rate_n = rate(w, s)

w, s = plot('tp-pdi', periodic, 'g', 1)
T, rate_t = rate(w, s)

ax = plt.subplot(gs[2])

ww = np.arange(0, 300, 0.1)
plt.plot(np.exp(-ww / 10) / (10 * (1 - np.exp(-300./10.))) , ww, 'c', label='T = 10K')
plt.plot(np.exp(-ww / 280) / (280 * (1 - np.exp(-300./280.))) , ww, 'm', label='T = 280K')
plt.fill_betweenx(ww, np.exp(-ww / 10) / (10 * (1 - np.exp(-300./10.))), alpha=0.5, color='c')
plt.fill_betweenx(ww, np.exp(-ww / 280) / (280 * (1 - np.exp(-300./280.))), alpha=0.5, color='m')
plt.yticks([])
plt.xticks([0., 0.01, 0.02])
plt.xlabel('$P_T(\\omega)$')
plt.xlim(0, 0.02)
plt.legend(handlelength=1.)

plt.subplots_adjust(wspace=0, hspace=0)

plt.savefig('bands.pdf', dpi=300)

plt.show()

#plt.plot(T, rate_n / rate_t)
#plt.show()
