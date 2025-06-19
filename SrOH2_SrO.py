from lib.SrOH2_volatile_models import *
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from lib.chemical_potentials import *


E_SrOH2 = chem_pot_SrOH2(1)
print("E SrOH2 in eV: ", E_SrOH2/ev2J_p_mol)

delta_E = E_SrO + E_DFT_H2O - E_SrOH2

print("delta_E in eV", delta_E/ev2J_p_mol)

T_lower_bound = 700 #K
T_upper_bound = 1200 #K

T_range = np.arange(T_lower_bound, T_upper_bound)
delta_G_list = []

for T in T_range:
    mu_H2O = chem_pot_H2O(T, E_DFT_H2O) + zpe_H2O
    delta_G = E_SrO + mu_H2O - chem_pot_SrOH2(T)
    delta_G_list.append(delta_G)

fig, ax = plt.subplots(layout="constrained")

delta_G_list = np.asarray(delta_G_list)
ax.plot(T_range, delta_G_list/ev2J_p_mol)


ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())


ax.set_xlim(T_lower_bound, T_upper_bound+1)

ax.set_xlabel("T [K]")
ax.set_ylabel("$\Delta_rG^*(T,p)$")


plt.show()
