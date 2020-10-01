# -*- coding: utf-8 -*-

import getpass
username = getpass.getuser()
from pathlib import Path

import numpy as np

# formatting for matplotlib
import matplotlib.pyplot as plt
plt.style.use(Path(f'C:\\Users\\{username}\\Dropbox (MIT)\\Research\\useful_code\\python\\paul_style.mplstyle'))
import matplotlib.ticker as ticker
import addcopyfighandler
plt.close(fig = 'all')


class species:
    """Describes each particle species in thermal equilibrium
    
    All energies and temperatures in eV.
    """
    
    def __init__(self, m_s, g_s, name_s, bf_val, is_QCD):
        self.m_s = m_s
        self.g_s = g_s
        self.name_s = name_s
        self.bf_val = bf_val # -1 if boson, +1 if fermion
        self.is_QCD = bool(is_QCD)
        self.lambda_QCD = 218e6
        
    def rho_s(self, T):
        """Energy density of species at given T"""
        if (self.m_s > T) or (self.is_QCD and T < self.lambda_QCD): #annihilated
            this_rho = 0
        elif self.bf_val == -1: #fermion
            this_rho = 7*self.g_s*np.pi**2*T**4/(8*30)
        elif self.bf_val == +1: #boson
            this_rho = self.g_s*np.pi**2*T**4/30
        return this_rho

class all_species:
    """A list of all species"""
    
    def __init__(self):
        # values of all species
        self.name_list = ["top quark", "bottom quark", "charm quark", "strange quark", "down quark", "up quark", "gluon", "taon", "muon", "electron", "tau neutrino", "muon neutrino", "electron neutrino", "W+", "W-", "Z0", "photon", "higgs boson"]
        self.mass_list = [173e9, 4e9, 1e9, 100e6, 5e6, 2e6, 0, 1777e6, 106e6, 511e3, 0.6, 0.6, 0.6, 80e9, 80e9, 91e9, 0, 125e9]
        self.g_list = [12, 12, 12, 12, 12, 12, 16, 4, 4, 4, 2, 2, 2, 3, 3, 3, 2, 1]
        self.bf_list = [-1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1]
        self.QCD_list = [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.n_species = len(self.name_list)
        
        # create list of all species
        self.species_list = []
        for s in range(self.n_species):
            this_species = species(self.mass_list[s], self.g_list[s], self.name_list[s], self.bf_list[s], self.QCD_list[s])
            self.species_list.append(this_species)
        
    def rho_sum(self, T):
        """Gets density of all species combined"""
        rho_out = 0
        for s in range(self.n_species):
            rho_out = rho_out + self.species_list[s].rho_s(T)
        return rho_out


# calculate rho and H vs T
all_spec = all_species()
n_pts = 1000
T_list = np.logspace(3, 12, n_pts)
rho_list = np.zeros(n_pts)
for T_ind in range(n_pts):
    rho_list[T_ind] = all_spec.rho_sum(T_list[T_ind])
H_list= np.sqrt(1.297e-25 * rho_list)


# plot rho/T^4 vs T
fig, ax = plt.subplots()
plt.loglog(T_list*1e-6, 30/np.pi**2*np.divide(rho_list, T_list**4), 'k-')
#plt.title('Problem 1(a)')
plt.xlabel('T [MeV]')
plt.ylabel('$\\rho/T^4$')
ax.xaxis.set_major_locator(ticker.LogLocator(numticks=4))
#ax.xaxis.set_minor_locator(ticker.LogLocator(numticks=10))
ax.yaxis.set_major_locator(ticker.LogLocator(numticks=4))
#ax.yaxis.set_minor_locator(ticker.LogLocator(numticks=15))


# plot H vs T
fig, ax = plt.subplots()
plt.loglog(T_list*1e-6, np.reciprocal(H_list), 'k-')
#plt.title('Problem 1(b)')
plt.xlabel('T [MeV]')
plt.ylabel('H [s]')
ax.xaxis.set_major_locator(ticker.LogLocator(numticks=4))
#ax.xaxis.set_minor_locator(ticker.LogLocator(numticks=10))
ax.yaxis.set_major_locator(ticker.LogLocator(numticks=4))
#ax.yaxis.set_minor_locator(ticker.LogLocator(numticks=15))