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


prefactor = 4*np.pi**3*const.G*(20)**0.5/(45*1.88e-29 *u.g*u.cm**-3)
#(2.725*u.K)**3
print(prefactor.decompose())