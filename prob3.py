# -*- coding: utf-8 -*-

import getpass
username = getpass.getuser()
from pathlib import Path
from astropy import constants as const
from astropy import units as u

import numpy as np
from scipy.optimize import root_scalar

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


def ionization_fract_eq(T):
    """calculates ionization fraction from temperature using Saha eq"""
    dE = 13.6
    m_e = 511e3
    factor = np.multiply(np.reciprocal(baryon_density(T)), np.multiply(np.power(m_e*T/2/np.pi, 3/2), np.exp(-np.divide(dE, T))))
    return 0.5*(np.sqrt(np.power(factor, 2) + 4 * factor) - factor)


def ion_z_eq_test(z):
    """returns ionization fraction at given redshift
    
    tests whether all the math is working by reproducing plot in textbook"""
    a = z_to_a(z)
    T = a_to_T(a)
    omega_b_h2 = 0.0222
    omega_m_h2 = 0.143
    return 1 - 123*ionization_fract_eq(T)
    
    
    

# define some constants
T_CMB = (2.725 * u.K * const.k_B).to(u.eV)
print(f"CMB temp today = {T_CMB}")
m_e = 511e3
dE = 13.6


# test plot
fig, ax = plt.subplots()
z_list = np.flip(np.logspace(4, 2, 1000))
X_e = ionization_fract_eq(a_to_T(z_to_a(z_list)))
X_e = X_e / X_e[-1]
plt.loglog(z_list, X_e, 'k-')
#plt.title('Problem 1(a)')
plt.xlabel('z')
plt.ylabel('$X_e$')
#ax.xaxis.set_major_locator(ticker.LogLocator(numticks=4))
#ax.xaxis.set_minor_locator(ticker.LogLocator(numticks=10))
#ax.yaxis.set_major_locator(ticker.LogLocator(numticks=4))
#ax.yaxis.set_minor_locator(ticker.LogLocator(numticks=15))
plt.ylim([1e-4,1.2])
plt.gca().invert_xaxis()


fig, ax = plt.subplots()
T = a_to_T(z_to_a(z_list))
plt.loglog(z_list, T, 'k-')
plt.gca().invert_xaxis()

fig, ax = plt.subplots()
y=np.power(m_e*T/2/np.pi, 3/2)
plt.loglog(z_list, y, 'k-')
plt.gca().invert_xaxis()

fig, ax = plt.subplots()
y = np.reciprocal(baryon_density(T))
plt.loglog(z_list, y, 'k-')
plt.gca().invert_xaxis()

fig, ax = plt.subplots()
y=np.exp(-np.divide(dE, T))
plt.loglog(z_list, y, 'k-')
plt.gca().invert_xaxis()