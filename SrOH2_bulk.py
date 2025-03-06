from lib.SrOH2_models import * 
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

x = 0.4 #molar fraction of Sr
x_O2 = 0.21
x_H2O = 0.08
P = 1 #atm
T_lower_bound = 700 #K
T_upper_bound = 1400 #K
T_range = np.arange(T_lower_bound, T_upper_bound)

p_O2 = x_O2 * P
p_H2O = x_H2O * P

fig, ax = plt.subplots(layout="constrained")
fig2, ax2 = plt.subplots(layout="constrained")

V_Sr, delta_G= case1(T_range)
ax.plot(T_range, V_Sr, label="case1")
ax2.plot(T_range, np.asarray(delta_G)/ev2J_p_mol, label="case1")

V_Sr, delta_G = case2(T_range)
ax.plot(T_range, V_Sr, label="case2")
ax2.plot(T_range, np.asarray(delta_G)/ev2J_p_mol, label="case2")

V_Sr, delta_G = case3(T_range)
ax.plot(T_range, V_Sr, label="case3")
ax2.plot(T_range, np.asarray(delta_G)/ev2J_p_mol, label="case3")

V_Sr, delta_G = case4(T_range)
ax.plot(T_range, V_Sr, label="case4")
ax2.plot(T_range, np.asarray(delta_G)/ev2J_p_mol, label="case4")


V_Sr, delta_G = case5(T_range)
ax.plot(T_range, V_Sr, label="case5")
ax2.plot(T_range, np.asarray(delta_G)/ev2J_p_mol, label="case5")

ax.set_title("Sr(OH)$_2$ formation")
ax.set_xlabel("T[K]")
ax.set_ylabel("[V\'\'\'$_{La}$]")


ax.set_xlim(left=T_lower_bound,right=T_upper_bound)
ax.set_ylim(0,)
ax.legend(loc="upper left")

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

ax2.set_title("Sr(OH)$_2$ formation")
ax2.set_xlabel("T[K]")
ax2.set_ylabel("$\Delta$G$^*$ [eV]")


ax2.set_xlim(left=T_lower_bound,right=T_upper_bound)
#ax2.set_ylim(0,)
ax2.legend(loc="upper left")
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())


ax2.legend()
plt.show()