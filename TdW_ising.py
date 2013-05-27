from numpy import array, exp, arange, roll, append
import matplotlib.pyplot as plt
from random import choice, randrange, random
from time import time

N = 20 #number of spins
J = 1. #coupling
n_therm = 100 #Thermalization of system
n = 400 
SPIN = [-1, 1]
#BETA = arange(0.1, 3.5, 0.1)
BETA = arange(0.1, 1, 0.1)
current_sample = []
magnetization = array([])
m_list = array([])
energy = array([])
e_list = array([])

def flip(spins, idx):
    """flips spin at position idx in spins"""
    flipped_spins = spins.copy()
    flipped_spins[idx] = -flipped_spins[idx]
    return flipped_spins

def iter_step(sample):
    """one step of metropolis algorithm, cf wikipedia Ising-model"""
    #Step 1: initialize random state
    spins = array([choice(SPIN) for i in range(N)])

    #Step 2: flip random spin
    idx = randrange(len(spins))
    flipped_spins = flip(spins, idx)

    #Step 3: calculate probability a
    if idx == len(spins)-1:
        H = -J * (spins[idx-1]*spins[idx] + spins[idx]*spins[0])
        H_flipped = -J * (flipped_spins[idx-1]*flipped_spins[idx] + 
                      flipped_spins[idx]*flipped_spins[0])
    else:
        H = -J * (spins[idx-1]*spins[idx] + spins[idx]*spins[idx+1])
        H_flipped = -J * (flipped_spins[idx-1]*flipped_spins[idx] + 
                          flipped_spins[idx]*flipped_spins[idx+1])

    a = exp(-beta*(H_flipped - H))
    if a < 0:
        a = 1.

    #Step 4: accept flipped_spins with prob. a
    if a >= random():
        sample = flipped_spins		
    else:
        sample = spins
    return sample

def update_m_e(sample):
    """updates magnetization, energy per sweep"""
    global magnetization, energy 
    magnetization = append(magnetization, array([sum(sample)]))
    energy = append(energy, -J*sum(sample*roll(sample, 1)))

def plot_magn_energy():
    """creates plot for energy and magnetization and saves it"""
    plt.plot(BETA, m_list, 'go', label="mean magnetization")
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
                update_m_e(sample)
        current_sample.append(sample) #add last sample to current array 
        m_list = append(m_list, magnetization.mean())
        e_list = append(e_list, energy.mean())
        t4 = time()
        print "time for beta = ", beta, " is: ", t4 - t3, "s"
    t2 = time()
    print "time taken: ", (t2 - t)/60, " min."
    plot_magn_energy()
    plt.show()
