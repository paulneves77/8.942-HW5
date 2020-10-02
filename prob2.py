# -*- coding: utf-8 -*-

import getpass
username = getpass.getuser()
from pathlib import Path
from astropy import constants as const
from astropy import units as u

import numpy as np
from scipy.special import zeta, gamma

# formatting for matplotlib
import matplotlib.pyplot as plt
plt.style.use(Path(f'C:\\Users\\{username}\\Dropbox (MIT)\\Research\\useful_code\\python\\paul_style.mplstyle'))
import matplotlib.ticker as ticker
import addcopyfighandler
plt.close(fig = 'all')


# calc number densities
photon_density = (zeta(3) * gamma(3) * (2.725 * u.K)**3 / np.pi**2)
photon_density = photon_density * const.k_B**3 * const.hbar**-3 * const.c**-3
print(f"photon density = {photon_density.to(u.cm**-3)}")
neutrino_density = 3 / 11* photon_density
print(f"neutrino density = {neutrino_density.to(u.cm**-3)}")


# calc temps
print(f"min temp = {(0.06 * u.eV / const.k_B).to(u.K)}")
print(f"max temp = {(0.12 * u.eV / const.k_B).to(u.K)}")


# neutrino temperature, a, z when became non-relativistic
T_nu_final = (4/11)**(1/3) * 2.725 * u.K
print(f"final neutrino temp = {T_nu_final}")

a_max = T_nu_final/(1390*u.K)
a_min = T_nu_final/(700*u.K)
print(f"max a = {a_max}")
print(f"min a = {a_min}")
z_max = 1/a_min-1
z_min = 1/a_max-1
print(f"max z = {z_max}")
print(f"min z = {z_min}")

print(f"avg z = {np.mean([z_max,z_min])}")
print(f"z unc = {np.ptp([z_max,z_min])}")


# calc min and max omega_nu*h^2
rho_crit = 3*(100 * u.km/u.s/u.Mpc)**2/(8*np.pi*const.G)
print(f"rho crit = {rho_crit.to(u.g*u.cm**-3)}")
prefactor = 112*u.cm**-3/(rho_crit*const.c**2)
print(f"inverse prefactor = {(prefactor**-1).to(u.eV)}")
print(f"neutrino fraction max = {(prefactor*0.06*u.eV).decompose()}")
print(f"neutrino fraction max = {(prefactor*0.12*u.eV).decompose()}")
