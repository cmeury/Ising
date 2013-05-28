from numpy import array, exp, arange, roll, append
import matplotlib.pyplot as plt
from random import choice, randrange, random
from time import time

#Parameters:
L = 20 #number of sites
d = 1 #dimension
N = L**d #number of spins
J = 1. #coupling
n_therm = 1000 #Thermalization of system
n = 4000 
BETA = arange(0.1, 3.5, 0.1)
#BETA = arange(0.1, 1, 0.1)
current_sample = array([])

magnetization = 0.
m_list = array([])
energy = 0.
e_list = array([])
N_magn = 0.
N_energy = 0.

class Spin:
    def __init__(self, L, d, J):
        """initializes a lattice of L**d spins with coupling constant J
           has spins, energy, magnetization"""
        SPIN = [-1, 1]
        self.spins = array([choice(SPIN) for i in range(L**d)])
        self.energy = -J * sum(self.spins * roll(self.spins, 1))
        self.magnetization = sum(self.spins)

    def flip(self, idx):
        """flips spin at position idx in spins"""
        flipped_spins = self.spins.copy()
        flipped_spins[idx] = -flipped_spins[idx]
        return flipped_spins

    def energy_diff(self, flipped_spins, idx):
        """returns difference of energy of flipped and unflipped spins"""
        N = len(self.spins)
        H = -J * (self.spins[idx-1]*self.spins[idx] + 
                  self.spins[idx]*self.spins[(idx+1)%N])
        H_flipped = -J * (flipped_spins[idx-1]*flipped_spins[idx] + 
                          flipped_spins[idx]*flipped_spins[(idx+1)%N])
        return  H_flipped - H

    def set_energy(self, J):
        self.energy = -J * sum(self.spins * roll(self.spins, 1))

    def set_magnetization(self):
        self.magnetization = sum(self.spins)

    def set_spins(self, new_spins):
        self.spins = new_spins
        self.set_energy(J)
        self.set_magnetization()

def iter_step(sample):
    """one step of metropolis algorithm, cf wikipedia Ising-model"""
    #Step 1: initialize random state
    spins = Spin(L, d, J) 

    #Step 2: flip random spin
    idx = randrange(len(spins.spins))
    flipped_spins = spins.flip(idx)

    #Step 3: calculate probability a
    energy_diff = spins.energy_diff(flipped_spins, idx)
    a = exp(-beta*(energy_diff))

    #Step 4: accept flipped_spins with prob. a
    if a >= random():
        return flipped_spins 
    else:
        return spins.spins

def update_energy_magnetization(sample):
    """updates magnetization, energy per sweep"""
    global magnetization, N_magn, energy, N_energy
    magnetization += abs(sum(sample))
    N_magn += 1.
    energy += -J*sum(sample*roll(sample, 1))
    N_energy += 1.

def plot_magn_energy():
    """creates plot for energy and magnetization and saves it"""
    plt.figure(1)

    plt.subplot(2,1,1)
    plt.plot(BETA, m_list, 'go', label="mean magnetization")
    plt.xlabel("beta")
    plt.ylabel("energy / magnetization")
    plt.legend(loc="center right")

    plt.subplot(2,1,2)
    plt.plot(BETA, e_list, 'ro', label="mean energy")
    plt.xlabel("beta")
    plt.ylabel("energy / magnetization")
    plt.legend(loc="center right")

    plt.savefig("Tdw_ising__mean_magn_and_energy.png")

if __name__ == "__main__":
    sample = 0
    t = time()
    for beta in BETA:
        t3 = time()
        #Thermalization process
        for i in range(n_therm):
            for k in range(N):
                sample = iter_step(sample) #Step 5
        #keep one sample per n sweeps
        for i in range(n):
            for k in range(N):
                sample = iter_step(sample) #Step 5
                update_energy_magnetization(sample)
        current_sample = append(current_sample, sample) 
        m_list = append(m_list, array([magnetization / N_magn]))
        e_list = append(e_list, energy.mean())
        t4 = time()
        print "time for beta = ", beta, " is: ", t4 - t3, "s"
    t2 = time()
    print "time taken: ", (t2 - t)/60, " min."
    plot_magn_energy()
    plt.show()
