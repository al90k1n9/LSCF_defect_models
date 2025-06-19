import matplotlib.pyplot as plt
from lib.SrO_humid_models import *
from matplotlib.ticker import  AutoMinorLocator

x = 0.4 #molar fraction of Sr
x_O2 = 0.21
x_H2O = 0.08
P = 1 #atm
T_lower_bound = 700
T_upper_bound = 800
T_range = np.arange(T_lower_bound,T_upper_bound,1) #K
#numpy imported chemical potentials, which is imported in humid models

p_O2 = x_O2 * P
p_H2O = x_H2O * P

fig,ax = plt.subplots(layout='constrained')
fig2, ax2 = plt.subplots(layout="constrained")


V_Sr, delta_G_range = case2(T_range, x, p_H2O, P=1)
ax.plot(T_range, V_Sr, label="[2H]$_{La}^{'}$")
ax2.plot(T_range, np.asarray(delta_G_range)/ev2J_p_mol, label="P=1 atm")


V_Sr, delta_G_range = case2(T_range, x, p_H2O, P=10)
ax2.plot(T_range, np.asarray(delta_G_range)/ev2J_p_mol, label="P=10 atm")

V_Sr, delta_G_range = case2(T_range, x, p_H2O, P=20)
ax2.plot(T_range, np.asarray(delta_G_range)/ev2J_p_mol, label="P=20 atm")

ax.axvline(x=158, color="black", ls="dashed")
ax.set_title("Humid conditions")
ax.set_xlabel("T[K]")
ax.set_ylabel("[V\'\'\'$_{La}$]")
ax.set_xlim(left=T_lower_bound,right=T_upper_bound)
ax.set_ylim(0,)
ax.legend(loc="upper right")
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())


ax2.axhline(y=0, color ="black", ls = "dashed")
ax2.axvline(x=158, color="black", ls="dashed")
ax2.set_title("")
ax2.set_xlabel("T[K]")
ax2.set_ylabel("$\Delta G^*$ [eV]")
ax2.set_xlim(left=T_lower_bound,right=T_upper_bound)
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())
ax2.legend()

bohr2m = 5.29177e-11
a_LSCF = 1.46415980775980e+01 * bohr2m
a_SrO = 3.14 * 1e-10 #m
sample_thickness = 20 * 1e-6 #m
specific_surface_area = 3.59*1e6 #m^2/m^3 <=> active surface area per unit volume of the electrode
volume_fraction_LSCF = 0.48

def yaxconvert(x):
    return x * volume_fraction_LSCF * a_SrO**2/(a_LSCF**3 * specific_surface_area)

def yaxinvert(x):
    return x/(volume_fraction_LSCF * a_SrO**2/(a_LSCF**3 * specific_surface_area))

secyax = ax.secondary_yaxis("right", functions=(yaxconvert, yaxinvert))
secyax.set_ylabel("blocked active surface area / total active surface area")


plt.show()
