import numpy as np
import sys
from plot import *

periodic = False

mol = sys.argv[1]

if mol == 'n' or mol == 'ne':
    labels = ['N-Phenyl-PDI Dynamics']
    molecule = 'N-Phenyl PDI'
elif mol == 't' or mol == 'te':
    labels = ['Tetraphenyl-PDI Dynamics']
    molecule = 'Tetraphenyl-PDI'
elif mol == 'b':
    labels = ['N-Phenyl-PDI Dynamics', 'Tetraphenyl-PDI Dynamics']
    molecule = 'PDI'
else:
    labels = sys.argv[2:]
    molecule = ''

plt.figure(figsize=(14,7))

if mol ==  'n':
    plt.xlim(0, 0.12)
    plt.ylim(0, 3.8)
if mol ==  't':
    plt.xlim(0, 0.3)
    plt.ylim(0, 3.5)
if mol == 'b':
    plt.xlim(0, 0.12)
    plt.ylim(0, 3.8)
if mol ==  'ne':
    plt.xlim(0, 0.3)
    plt.ylim(0, 4.2)
if mol == 'te':
    plt.xlim(0, 0.8)
    plt.ylim(0, 3.6)
if mol == 'be':
    plt.xlim(0, 0.3)
    plt.ylim(0, 4.2)

for k, filename in enumerate(sys.argv[2:]):

    times, pops = load(filename)
    
    assert pops.shape[1] in [i**2 for i in range(150)]
    
    term_0 = np.zeros_like(times)
    term_1 = np.zeros_like(times)
    sep = np.zeros_like(times)
    
    nsites = int(np.sqrt(pops.shape[1]))
    origin = np.where(pops[0])[0][0] / (nsites + 1)
    
    for t in range(len(times)):
    
        for i in range(nsites):
    
            for j in range(nsites):
      
                if periodic: 
                    d = min(abs(i - j), nsites - abs(i - j))
                else:
                    d = abs(i - j) 
    
                sep[t] += d * pops[t, nsites * i + j]
    
                if (i + d) % nsites == j:
                    left = float(i)
                elif (j + d) % nsites == i:
                    left = float(j)
                else:
                    raise ValueError
    
                center = ((left + d / 2.) % nsites)
    
                if periodic:
                    rms = min((center - origin)**2, (nsites - abs(center - origin))**2)
                else:
                    rms = abs(center - origin)**2
    
                #com[t] += rms * pops[t, nsites * i + j]
                term_0[t] += pops[t, nsites * i + j] * (0.5 * (i + j))**2 
                term_1[t] += pops[t, nsites * i + j] * 0.5 * (i + j)
                
    com = term_0 - term_1**2
    com = np.sqrt(com)

    plt.plot(times, com, label=labels[k])

if mol == 'n' or mol == 'b':
    plt.plot(times, 0.336687 + 24.765 * times, label='Linear approximation to N-Phenyl PDI dynamics')
if mol == 't' or mol == 'b':
    plt.plot(times, 0.47 + 8.4 * times, label='Linear approximation to Tetraphenyl PDI dynamics')
plt.xlabel('Time (ps)')
plt.ylabel('$\\langle\\Delta x(t)^2\\rangle$')
plt.legend(loc='upper left')
plt.title('RMS distance of center of mass from exciton origin in units of intermolecular spacings, '+molecule)
plt.savefig('fig.pdf', dpi=300, bbox='tight')
plt.show()

#plt.plot(times, sep)
#plt.xlabel('Time (ps)')
#plt.title('Electron-hole separation, in units of intermolecular spacings, N-Phenyl-PDI')
#plt.show()
