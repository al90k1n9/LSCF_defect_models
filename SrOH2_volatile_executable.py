from lib.SrOH2_volatile_models import *
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator



x = 0.4 #molar fraction of Sr
x_O2 = 0.21
x_H2O = 0.08
P = 1 #atm


data = np.genfromtxt("./lib/sroh2_factsage.csv", delimiter=";")
T_range = data[100::100,0]
print(np.shape(T_range))
V_Sr, delta_G= case1(T_range)

fig, ax = plt.subplots(layout="constrained")

ax.plot(T_range, V_Sr)

ax.set_xlabel("T [K]")
ax.set_ylabel("[V$_{Sr}$\'\'\']")
ax.set_xlim(T_range[0], T_range[-1])
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())



fig5, ax5 = plt.subplots(layout="constrained")
ax5.plot(T_range, delta_G/ev2J_p_mol)

ax5.set_xlabel("T [K]")
ax5.set_ylabel("$\Delta$G$^*$ [eV]")

ax5.set_xlim(T_range[0], T_range[-1])

ax5.xaxis.set_minor_locator(AutoMinorLocator())
ax5.yaxis.set_minor_locator(AutoMinorLocator())


fig.savefig("sroh2_volatile_case.png", format="png", dpi=300, transparent=True)
fig5.savefig("sroh2_volatile_delta_G.png", format="png", dpi=300, transparent=True)

plt.show()
