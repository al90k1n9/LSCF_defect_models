import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from lib.SrO_hydroxylated_models import *
from lib.auxilliary_functions import *

x0 = 0.4 #molar fraction of Sr
x_O2 = 0.21
x_H2O = 0.08
P = 1 #atm
T_lower_bound = 400
T_upper_bound = 1299
T_range = np.arange(T_lower_bound,T_upper_bound,1) #K
#numpy imported chemical potentials, which is imported in humid models

p_O2 = x_O2 * P
p_H2O = x_H2O * P


fig,ax = plt.subplots(layout='constrained')
axinset = ax.inset_axes([0.45,0.15,0.5,0.5])
#axinset.set_facecolor("none")
fig2, ax2 = plt.subplots(layout="constrained")

#fig3, ax3 = plt.subplots(layout="constrained")
#
#fig4, ax4  = plt.subplots(layout="constrained")


#V_Sr, delta_G, theta_list = case1(T_range, x0, p_O2, p_H2O, P)
#ax.plot(T_range, V_Sr, label="case 3.3")
#axinset.plot(T_range, V_Sr, label="case 3.3")
#ax2.plot(T_range, delta_G/ev2J_p_mol, label="case 3.3")
#ax4.plot(T_range, theta_list, label="case1")

V_Sr, delta_G, theta_list = case2(T_range, x0, p_H2O, P)
ax.plot(T_range, V_Sr, label="case 3.3")
axinset.plot(T_range, V_Sr/x0)
axinset.plot(T_range, theta_list, label="$\\theta_{H_2O}$")
axinset.legend(facecolor="none")
axinset.legend(facecolor="none")
ax2.plot(T_range, delta_G/ev2J_p_mol, label="case 3.3")
print("case 2 ", delta_G[1]/ev2J_p_mol)

#ax3.plot(T_range, np.asarray(V_Sr)/0.4*100, label="SrO")
#ax4.plot(T_range, theta_list, label="case1")

#V_Sr, delta_G, theta_list = case3(T_range, x0, p_O2, p_H2O, P)
#ax.plot(T_range, V_Sr, label="case 3")
#ax2.plot(T_range, delta_G/ev2J_p_mol, label="case3")
#ax4.plot(T_range, theta_list, label="case1")


#===================================================
#PLOTTING OPTIONS
ax.set_xlabel("T[K]")
ax.set_ylabel("$[V_{La}]_{eq}^{\'\'\'} = x_{eq}$")


ax.set_xlim(left=T_lower_bound,right=T_upper_bound+1)
ax.set_ylim(0,)
axinset.set_xlim(left=T_lower_bound,right=T_upper_bound+1)
axinset.set_ylim(-0.1,1.1)
ax.legend(loc="upper right", facecolor="none")

bohr2m = 5.29177e-11
a_LSCF = 1.46415980775980e+01 * bohr2m
a_SrO = 3.14 * 1e-10 #m
sample_thickness = 20 * 1e-6 #m
specific_surface_area = 3.59*1e6 #m^2/m^3 <=> active surface area per unit volume of the electrode
volume_fraction_LSCF = 0.48

#def yaxconvert(x):
#    return x * volume_fraction_LSCF * a_SrO**2/(a_LSCF**3 * specific_surface_area)
#
#def yaxinvert(x):
#    return x/(volume_fraction_LSCF * a_SrO**2/(a_LSCF**3 * specific_surface_area))
def yaxconvert(x):
    return x * 100/0.4

def yaxinvert(x):
    return x *0.4/100



secyax = ax.secondary_yaxis("right", functions=(yaxconvert, yaxinvert))
secyax.yaxis.set_minor_locator(AutoMinorLocator())
secyax.set_ylabel("% of total Sr content $\\frac{100\\cdot x_{eq}}{x_0}$")

#secyaxinset = axinset.secondary_yaxis("right", functions=(yaxconvert, yaxinvert))
#secyaxinset.yaxis.set_minor_locator(AutoMinorLocator())

ax2.set_xlabel("T[K]")
ax2.set_ylabel("$\Delta_rG^*(T,p)$ [eV]")

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
axinset.xaxis.set_minor_locator(AutoMinorLocator())
axinset.yaxis.set_minor_locator(AutoMinorLocator())
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())


ax2.set_xlim(left=T_lower_bound,right=T_upper_bound)
#ax2.set_ylim(0,)
ax2.legend(loc="upper right", facecolor="none")

#T_range = np.arange(T_lower_bound, T_upper_bound+200, 1)
#theta_list = surface_coverage_H2O(T_range, x_H2O=0.08, E_ads= E_ads, P=1)
#ax3.plot(T_range, theta_list*100, color="black", label="surface coverage")
#
#ax3.set_xlabel("T [K]")
#ax3.set_ylabel("% of total Sr")
#
#ax3.xaxis.set_minor_locator(AutoMinorLocator())
#ax3.yaxis.set_minor_locator(AutoMinorLocator())
#
#ax3.set_xlim(T_lower_bound, T_upper_bound)
#ax3.legend()
#
#def percent2surf_coverage(x):
#    return x/100
#
#def surf_coverage2percent(x):
#    return x*100
#
#ax32 = ax3.secondary_yaxis('right', functions=(percent2surf_coverage, surf_coverage2percent))
#ax32.set_ylabel("$\\theta$")
#ax32.yaxis.set_minor_locator(AutoMinorLocator())
#
##fig3.savefig("hydroxylated_case.png", dpi=300, transparent=True)
#
#ax4.set_xlabel("T [K]")
#ax4.set_ylabel("$\\theta$")
#
#ax4.xaxis.set_minor_locator(AutoMinorLocator())
#ax4.yaxis.set_minor_locator(AutoMinorLocator())
#
#ax4.set_xlim(T_lower_bound, T_upper_bound)
##ax4.set_ylim(0,)

fig.savefig("figs/hydroxylated_cond_delta_G.svg", format="svg", dpi=300, transparent=True)
fig2.savefig("figs/hydroxylated_cond_VSr.svg", format="svg", dpi=300, transparent=True)

plt.show()
