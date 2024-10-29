import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

ev2J = 1.60219*10**(-19)
lattice_SrO = 3.8 #Angstrom
gamma_SrO = 0.428 #eV per conventional surface

area_SrO = lattice_SrO**2 * 10**-20 #m2 
gamma_SrO *= area_SrO #eV/m2

print(gamma_SrO*ev2J)

W = np.arange(0.01, 5, 0.01)

interface_energy = (area_SrO * (2*gamma_SrO - W)/ev2J)


fig, ax = plt.subplots(layout="constrained")
ax.plot(W, interface_energy)
ax.set_xlabel("W [J/m2]")
ax.set_ylabel("interface/surface energy correction [eV]")
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlim(0,5)

plt.show()