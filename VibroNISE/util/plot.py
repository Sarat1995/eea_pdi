import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.style
matplotlib.style.use('classic')

def load(filename):

    with open(filename, 'r') as f:
    
        lines = f.readlines()
    
    data = []
    
    for line in lines:
    
        line = line.rstrip()
        line = line.rstrip(',')
    
        line = line.split(', ')
    
        line = [float(x) for x in line]
    
        data.append(line)
    
    data = np.array(data)
    
    times = data[:, 0]
    pops = data[:, 1:]

    return times, pops

if __name__ == "__main__":

    filename = sys.argv[1]

    times, pops = load(filename)

    pops = np.concatenate((pops, pops[:, [0]]), axis=1)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    Yn, Xn = np.meshgrid(times, np.arange(pops.shape[1]))

    density = ax.pcolormesh(Yn, Xn, pops.transpose(), vmin=0, vmax=1.)
    fig.colorbar(density, ax=ax)

    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(times, np.sum(pops[:, :-1], axis=1))
    plt.show()
