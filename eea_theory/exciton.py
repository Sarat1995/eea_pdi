"""exciton.py"""

import numpy as np
import os
from scipy.special import binom
import matplotlib.pyplot as plt

from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

class Exciton:

    def __init__(self, molecule='n-ph-pdi', beta=None, periodic=False):

        self.molecule = molecule
        assert self.molecule in ['n-ph-pdi', 'tp-pdi']

        self.periodic = periodic

        self.w = None
        self.v = None
        self.load_wavefunction(beta)
        #self.calculate_basis_positions()

#    def calculate_basis_positions(self):
#
#        #one_particle = self.nsites * self.nvib
#        #two_particle = int(binom(self.nsites, 2)) * self.nvib * (self.nvib - 1)
#
#        #self.nbasis = one_particle + two_particle
#
#        self.pos = []
#
#        # Calculate positions of one-particle basis states.
#        for i in range(self.nsites):
#            self.pos += [i] * self.nvib
#
#        # Calculate positions of two-particle basis states.
#        #for i in range(self.nsites):
#        #    self.pos += [i] * (int(binom(self.nsites, 2)) * self.nvib * (self.nvib - 1) / self.nsites)

    def wavefunction(self):

        return self.loaded_wavefunction()

    def load_wavefunction(self, beta):

        if self.periodic:
            p = '_periodic'
        else:
            p = ''
        
        if os.path.exists('../VibroNISE/runs/'):
            dirname = '../VibroNISE/runs/'
            dirname = os.path.join(dirname, self.molecule + p) 
        else:
            dirname = ''

        try:
            if rank == 0:
                self.w = np.load(self.molecule + p + '_evals.npy')
        except:
            if rank == 0:
                self.w = np.genfromtxt(os.path.join(dirname, 'SEigenValues.dat'), delimiter=',')[:, :-1].ravel()
                np.save(self.molecule + p + '_evals.npy', self.w)

        if rank == 0:
            self.w -= np.min(self.w)

        self.w = comm.bcast(self.w, root=0)

        # Boltzmann factors.
        boltzmann = np.exp(-beta * self.w)
  
        # Canonical partition function.
        Z =  np.sum(boltzmann)
        boltzmann /= Z

        inds = np.where(boltzmann >= 1.e-4)

        self.w = self.w[inds[0]]

        try:
            if rank == 0:
                self.v = np.load(self.molecule + p + '.npy')
        except:
            if rank == 0:
                self.v = np.genfromtxt(os.path.join(dirname, 'SEigenVectors.dat'), delimiter=',')[:, :-1].transpose()
                self.v = self.v[:, inds[0]]
                np.save(self.molecule + p + '.npy', self.v)

        if rank == 0:
            self.v = self.v[:, inds[0]]
            
        self.v = comm.bcast(self.v, root=0)

        try:
            with open(os.path.join(dirname, 'InputDirect1D.dat'), 'r') as f:
                lines = f.readlines()
        except:
            with open('InputDirect1D_' + self.molecule + p + '.dat', 'r') as f:
                lines = f.readlines()

        self.nvib = int(lines[lines.index('VMax\n') + 1]) + 1
        self.nsites = int(lines[lines.index('UnitTotal\n') + 1])

    def loaded_wavefunction(self):

        nstates = self.nvib * self.nsites

        #print "Norm:"
        #print np.linalg.norm(self.v, axis=0)

        v = np.copy(self.v[:nstates])

        for i in range(1, self.nvib):
            v[i::self.nvib] = 0.

        #print "New norm:"
        #print np.linalg.norm(v, axis=0)

        with np.errstate(invalid='ignore'):
            v /= np.linalg.norm(v, axis=0)

        #for i in range(self.nsites):
        #    #print (v[::self.nvib, i])
        #    print np.sum(self.w[i])
        #    print np.sum(v[::self.nvib, i])
        #    plt.plot(v[::self.nvib, i])
        #    plt.show()

        return v

    def sin_wavefunction(self, N, aggregate):

        raise NotImplementedError

        # Define domain.
        x = np.linspace(0, 1, num=N+2)[1:-1]

        # Calculate particle in a box wavefunction.
        if aggregate == 'j':
            psi = np.sin(np.pi * x)
        else:
            psi = np.sin(np.pi * x * N)

        # Normalize wavefunction.
        psi /= np.linalg.norm(psi)

        return psi

if __name__ == "__main__":

    exciton = Exciton()
    print(exciton.wavefunction().shape)
