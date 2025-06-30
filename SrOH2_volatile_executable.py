"""Models for the creation of SrOH2(g) volatile species.
"""

from lib.SrOH2_volatile_models import *
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os

local_path = os.path.dirname(os.path.abspath(__file__))
local_path += "/"


x = 0.4 #molar fraction of Sr
x_O2 = 0.21
x_H2O = 0.08
P = 1 #atm


T_lower_bound = 700
T_upper_bound = 1100
T_range = np.arange(T_lower_bound, T_upper_bound)



delta_oxygen_parameters = delta_oxygen_interpolater(plot=0)
V_Sr, delta_G, delta_oxygen_list, Ni_H2O_list = case1(T_range, delta_oxygen_parameters=delta_oxygen_parameters)

fig, ax = plt.subplots(layout="constrained")

ax.plot(T_range, V_Sr)

ax.set_xlabel("T [K]")
ax.set_ylabel("$[V_{La}\'\'\']_{eq} = x_{eq}$")
ax.set_xlim(T_range[0], T_range[-1])
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

def yaxconvert(x):
    return x * 100/0.4

def yaxinvert(x):
    return x *0.4/100

secyax = ax.secondary_yaxis("right", functions=(yaxconvert, yaxinvert))
secyax.set_ylabel("% of initial Sr content $\\frac{100 \cdot x_{eq}}{x_0}$")
secyax.yaxis.set_minor_locator(AutoMinorLocator())


fig2, ax2 = plt.subplots(layout="constrained")
ax2.plot(T_range, delta_G/ev2J_p_mol)

ax2.set_xlabel("T [K]")
ax2.set_ylabel("$\Delta_r$G$^*(T,p)$ [eV]")

ax2.set_xlim(T_range[0], T_range[-1])

ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())

fig3, ax3 = plt.subplots(layout="constrained")
ax3.plot(T_range, Ni_H2O_list)

ax3.set_xlabel("T [K]")
ax3.set_ylabel("initial # H_2O molecules/LSCF unit cell")
ax3.set_xlim(T_range[0], T_range[-1]+1)
ax3.xaxis.set_minor_locator(AutoMinorLocator())
ax3.yaxis.set_minor_locator(AutoMinorLocator())


fig4, ax4 = plt.subplots(layout="constrained") #for sensitivity analysis
Ni_H2O, solutions, x_H2O_list = Ni_H2O_sensitivity()
ax4.plot(x_H2O_list, solutions)

ax4.set_xlabel("$x_{H_2O}$")
ax4.set_ylabel("$[V_{La}\'\'\']_{eq} = x_{eq}$")
#ax4.set_xlim()
ax4.set_ylim(0,)
ax4.yaxis.set_minor_locator(AutoMinorLocator())



ax4.set_xscale("log")


fig.savefig(local_path + "figs/sroh2_volatile_case.png", format="png", dpi=300, transparent=True)
fig2.savefig(local_path + "figs/sroh2_volatile_delta_G.png", format="png", dpi=300, transparent=True)

plt.show()
