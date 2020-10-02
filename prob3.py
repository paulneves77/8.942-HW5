# -*- coding: utf-8 -*-

import getpass
username = getpass.getuser()
from pathlib import Path
from astropy import constants as const
from astropy import units as u

import numpy as np
from scipy.optimize import root_scalar
from functools import partial

# formatting for matplotlib
import matplotlib.pyplot as plt
plt.style.use(Path(f'C:\\Users\\{username}\\Dropbox (MIT)\\Research\\useful_code\\python\\paul_style.mplstyle'))
import matplotlib.ticker as ticker
import addcopyfighandler
plt.close(fig = 'all')


def z_to_a(z):
    """converts expansion parameter value to redshift value"""
    return np.divide(1, z+1)


def a_to_T(a):
    """gets the temperature as a function of expansion parameter"""
    return np.divide(0.00023482, a)


def baryon_density(T):
    """gets baryon number density as a function of temperature"""
    return 1e-9 * np.power(T, 3)


def T_to_Xe(T):
    """calculates ionization fraction from temperature using Saha eq"""
    dE = 13.6
    m_e = 511e3
    factor = np.multiply(np.reciprocal(baryon_density(T)), np.multiply(np.power(m_e*T/2/np.pi, 3/2), np.exp(-np.divide(dE, T))))
    return 0.5*(np.sqrt(np.power(factor, 2) + 4 * factor) - factor)


def ion_z_eq_full(obh2, z):
    """returns ionization fraction at given redshift
    
    tests whether all the math is working by reproducing plot in textbook"""
    a = z_to_a(z)
    T = a_to_T(a)
    omh2 = 0.143
    return 1-123*T_to_Xe(T)*(obh2/0.022)*(0.14/omh2)**(1/2)*((1+z)/1000)**(3/2)*(1+(1+z)*0.14/(3360*omh2))**(-1/2)


# define some constants
T_CMB = (2.725 * u.K * const.k_B).to(u.eV)
print(f"CMB temp today = {T_CMB}")
m_e = 511e3
dE = 13.6


# test plot
fig, ax = plt.subplots()
z_list = np.flip(np.logspace(4, 2, 1000))
X_e = T_to_Xe(a_to_T(z_to_a(z_list)))
#X_e = X_e / X_e[-1]
plt.loglog(z_list, X_e, 'k-')
#plt.title('Problem 1(a)')
plt.xlabel('z')
plt.ylabel('$X_e$')
plt.ylim([1e-4,1.2])
plt.gca().invert_xaxis()


# plot z recombine vs omega_b*h^2
obh2_list = np.linspace(0.01, 0.05, 100)
z_recombine = np.zeros(len(obh2_list))
for ind in range(len(obh2_list)):
    obh2 = obh2_list[ind]
    ion_z_eq = partial(ion_z_eq_full, obh2)
    this_root = root_scalar(ion_z_eq, x0=1110, x1=1e3)
    z_recombine[ind] = this_root.root
    
fig, ax = plt.subplots()
plt.plot(obh2_list, z_recombine, 'k-')
#plt.title('Problem 1(a)')
plt.xlabel('$\\Omega_b h^2$')
plt.ylabel('$z_{recombine}$')
#plt.xlim([1e-2,5e-2])