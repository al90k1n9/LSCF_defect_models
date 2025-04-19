import matplotlib.pyplot as plt
from lib.SrO_humid_models import *
from matplotlib.ticker import  AutoMinorLocator

x0 = 0.4 #molar fraction of Sr
x_O2 = 0.21
x_H2O = 0.08
P = 1 #atm
T_lower_bound = 700
T_upper_bound = 1299
T_range = np.arange(T_lower_bound,T_upper_bound,1) #K
#numpy imported chemical potentials, which is imported in humid models


fig,ax = plt.subplots(layout='constrained')
fig2, ax2 = plt.subplots(layout="constrained")

axinset = ax.inset_axes([0.4, 0.1, 0.5, 0.5])

V_Sr, delta_G_range = case1(T_range, x0, x_O2, x_H2O, P=1)
#keep in mind that the pH2 is determined inside the case1 function by pO2 and pH2O by considering the equilibrium between these three gases.
ax.plot(T_range, V_Sr, label="H$_{2(g)}$")
axinset.plot(T_range, V_Sr, label="H$_{2(g)}$")
ax2.plot(T_range, np.asarray(delta_G_range)/ev2J_p_mol, label="H$_{2(g)}$")

V_Sr, delta_G_range = case2(T_range, x0, x_H2O, P=1)
ax.plot(T_range, V_Sr, label="(2H)$_{La}^{'}$")
ax2.plot(T_range, np.asarray(delta_G_range)/ev2J_p_mol, label="(2H)$_{La}^{'}$")


#V_Sr, delta_G_range = case3(T_range, x0=0.4, x_O2 = 0.21, x_H2O = 0.08, P=1)
#ax.plot(T_range, V_Sr, label="$\\frac{1}{2}$H$_{2(g)}$+[H]$_{La}^{''}$")
#ax2.plot(T_range, np.asarray(delta_G_range)/ev2J_p_mol, label="case 3")
#
#
#V_Sr, delta_G_range = case4(T_range, x0, x_O2, x_H2O, P=1)
#ax.plot(T_range, V_Sr, label="H$_{2(g)}$, Sr$_{bulk}$")
#ax2.plot(T_range, np.asarray(delta_G_range)/ev2J_p_mol, label="case 4")

V_Sr, delta_G = case5(T_range)
#ax.plot(T_range, V_Sr, label="case 4")
#ax2.plot(T_range, delta_G/ev2J_p_mol, label="case 5")


fig3, ax3 = plt.subplots(layout = "constrained")


def func_case5(x, T, delta_G, volume_fraction_gas = 0.5):
    volume_gas = volume_fraction_gas/(1-volume_fraction_gas) * (acell_LSCF_slab/2)**3
    numerator = (x0 + 2*x)**2 * x *x*kB*T/(P*volume_gas)
    denominator = (1-x0-2*x)**2 * (x0-x)
    return np.exp(-delta_G/(R*T))-numerator/denominator

index = 200
f = lambda x:func_case5(x, T_range[index], delta_G[index], volume_fraction_gas=0.5)

xlist = np.linspace(0, 1e-14, 100)
ylist = []
for x in xlist:
    ylist.append(f(x))
ax3.plot(xlist, ylist)
ax3.axhline(y=0, color="black", ls="dashed")

#===================================================
#PLOTTING OPTIONS
ax.set_xlabel("T[K]")
ax.set_ylabel("$x_{eq} = [V\'\'\'_{La}$]")
ax.set_xlim(left=T_lower_bound,right=T_upper_bound+1)
ax.set_ylim(0,)
ax.legend(loc="upper right", facecolor="none")
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

def yaxconvert(x):
    return x * volume_fraction_LSCF * a_SrO**2/(a_LSCF**3 * specific_surface_area)

def yaxinvert(x):
    return x/(volume_fraction_LSCF * a_SrO**2/(a_LSCF**3 * specific_surface_area))

axinset.set_xlim(left=T_lower_bound,right=T_upper_bound+1)
axinset.set_ylim(1.8e-8, 6e-8)
axinset.xaxis.set_minor_locator(AutoMinorLocator())
axinset.yaxis.set_minor_locator(AutoMinorLocator())
axinset.set_facecolor("none")

secyax_inset = axinset.secondary_yaxis("right", functions=(yaxconvert, yaxinvert))
secyax_inset.yaxis.set_minor_locator(AutoMinorLocator())



ax2.set_title("")
ax2.set_xlabel("T[K]")
ax2.set_ylabel("$\Delta G^*(T, p)$ [eV]")
ax2.set_xlim(left=T_lower_bound,right=T_upper_bound+1)
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())
ax2.legend(facecolor="none")


bohr2m = 5.29177e-11
a_LSCF = 1.46415980775980e+01 * bohr2m
a_SrO = 3.14 * 1e-10 #m
sample_thickness = 20 * 1e-6 #m
specific_surface_area = 3.59*1e6 #m^2/m^3 <=> active surface area per unit volume of the electrode
volume_fraction_LSCF = 0.48

secyax = ax.secondary_yaxis("right", functions=(yaxconvert, yaxinvert))
secyax.set_ylabel("% of initial Sr content $\\frac{100 \cdot x_{eq}}{x_0}$")
secyax.yaxis.set_minor_locator(AutoMinorLocator())

fig.savefig("figs/humid_conditions_VSr.svg", format="svg", dpi=300, transparent=True)
fig2.savefig("figs/humid_conditions_delta_G.svg", format="svg", dpi=300, transparent=True)


plt.show()
