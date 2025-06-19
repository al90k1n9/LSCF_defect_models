"""SrO formation in dry air conditions

Different cases can be run by passing only the temperature window as arguments, with the rest being the default values, to get the results in the manuscript.
"""

import matplotlib.pyplot as plt 

from lib.SrO_dry_air_models import *
#chemical potentials are imported in models
#numpy imported chemical potentials
from matplotlib.ticker import  AutoMinorLocator
import os

local_path = os.path.dirname(os.path.abspath(__file__))
local_path += "/"

x = 0.4 #molar fraction of Sr
x_O2 = 0.21
P = 1 #atm
T_lower_bound = 600 
T_upper_bound = 1200
T_range = np.arange(T_lower_bound,T_upper_bound,1) #K

p_O2 = x_O2 * P


default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

fig,ax = plt.subplots(layout='constrained')

axinset1 = ax.inset_axes([0.07,0.43,0.3,0.3])
axinset2 = ax.inset_axes([0.65,0.43,0.3,0.3])
axinset1.set_facecolor("none")
axinset2.set_facecolor("none")

fig2, ax2 = plt.subplots(layout ="constrained")

fig3, ax3 = plt.subplots(layout="constrained")

fig4, ax4 = plt.subplots(layout="constrained")

fig5, ax5 = plt.subplots(layout="constrained")

fig6, ax6 = plt.subplots(layout="constrained")

ylist, delta_G_list = case1(T_range, x, p_O2, P)
plt_element_case1 = ax.plot(T_range, ylist, label ="case 1")
axinset1.plot(T_range, ylist)
plt_element_case1_dG = ax2.plot(T_range, delta_G_list/ev2J_p_mol, label="case 1")

index = np.where(T_range==973)[0][0]
K_eq = (ylist * (0.4 + 2*ylist)**2)/((0.4-ylist)*(1-0.4-ylist)**2 * np.sqrt(0.21))
print("case1  (index, temp [K], delta_G [eV], K(delta_G),  solution [molar concentration], K (solution x), diff K")
print(index, T_range[index], delta_G_list[index]/ev2J_p_mol, np.exp(-delta_G_list[index]/(R*T_range[index])), ylist[index]/0.4, K_eq[index], K_eq[index]-np.exp(-delta_G_list[index]/(R*T_range[index])))
ax5.plot(T_range, K_eq-np.exp(-delta_G_list/(R*T_range)), label="case1")


ylist, delta_G_list = case2(T_range, x, p_O2, P)
plt_element_case2, = ax.plot(T_range, ylist, label = "case 2")
axinset2.plot(T_range, ylist, color = default_colors[1])
plt_element_case2_dG = ax2.plot(T_range, np.asarray(delta_G_list)/ev2J_p_mol, label="case 2")

delta_oxygen_parameters = delta_oxygen_interpolater(plot=0)
#delta_oxygen_parameters = np.asarray([0,0])
delta_oxygen_inversion_temperature = -delta_oxygen_parameters[1]/delta_oxygen_parameters[0]


ylist, delta_G_list, delta_oxygen_list = case3(T_range, x=0.4, delta_oxygen_parameters=delta_oxygen_parameters)
plt_element_case3, =ax.plot(T_range, ylist, label ="case 3")
ax.axvline(x=delta_oxygen_inversion_temperature, color="black", ls="dashed")
plt_element_case3_dG = ax2.plot(T_range, np.asarray(delta_G_list)/ev2J_p_mol, label="case 3")
ax6.plot(T_range, delta_oxygen_list, label="case 3", color= default_colors[2])

index = np.where(T_range==973)[0][0]

delta_oxygen = delta_oxygen_parameters[0]*T_range[index] + delta_oxygen_parameters[1]
print(delta_oxygen, 0.00792)
K_eq = ylist * (ylist + delta_oxygen_list) /(((3-delta_oxygen_list)-ylist)*(0.4-ylist))
print("case3  (index, temp [K], delta_G [eV], K(delta_G),  solution [molar concentration], K (solution x), diff K ")
print(index, T_range[index], delta_G_list[index]/ev2J_p_mol, np.exp(-delta_G_list[index]/(R*T_range[index])), ylist[index]/0.4, K_eq[index], K_eq[index]-np.exp(-delta_G_list[index]/(R*T_range[index])))
ax5.plot(T_range, K_eq-np.exp(-delta_G_list/(R*T_range)), label="case3")

ylist, delta_G_list, delta_oxygen_list = case4(T_range, x = 0.4, delta_oxygen_parameters=delta_oxygen_parameters)
plt_element_case4, = ax.plot(T_range, ylist, label ="case 4")
plt_element_case4_dG = ax2.plot(T_range, np.asarray(delta_G_list)/ev2J_p_mol, label="case 4")
ax6.plot(T_range, delta_oxygen_list, label="case 6", color=default_colors[3])


#ylist, delta_G_list = case5(T_range, x, p_O2, P)
#plt_element_case5 = ax.plot(T_range, ylist, label ="case 5")
##axinset2.plot(T_range, ylist, color=default_colors[4])
#plt_element_case5_dG = ax2.plot(T_range, np.asarray(delta_G_list)/ev2J_p_mol, label="case 5")
#
#ylist, delta_G_list = case6(T_range, x=0.4)
#plt_element_case6 = ax.plot(T_range, ylist, label ="case 6")
#plt_element_case6_dG = ax2.plot(T_range, np.asarray(delta_G_list)/ev2J_p_mol, label="case 6")

x_O2_range = np.linspace(0.1, 1, 10)
for x_O2 in x_O2_range:
    ylist, delta_G_list = case1(T_range, x_O2= x_O2)
    ax3.plot(T_range, ylist, label="x$_{O_2}$ = " + str(round(x_O2, 2)))

total_p_range = np.linspace(3, 30, 10)
for total_p in total_p_range:
    ylist, delta_G_list = case1(T_range, P=total_p)
    ax4.plot(T_range, ylist, label="P=" + str(int(total_p))+" atm")


ax.set_xlabel("T[K]")
ax.set_ylabel("$[V\'\'\'_{La}]_{eq}=x_{eq}$")


ax2.set_xlabel("T[K]")
ax2.set_ylabel("${\Delta}G^*(T,p)$ [eV]")

ax3.set_xlabel("T[K]")
ax3.set_ylabel("$[V\'\'\'_{La}]_{eq}=x_{eq}$")

ax4.set_xlabel("T[K]")
ax4.set_ylabel("$[V\'\'\'_{La}]_{eq}=x_{eq}$")

ax5.axvline(x=973, color="black", ls="dashed")
ax5.set_xlabel("T [K]")
ax5.set_ylabel("delta G - RT ln K ")
ax5.legend()

ax.set_xlim(left=T_lower_bound,right=T_upper_bound+1)
ax.set_ylim(0,)
ax.legend(loc="upper left", facecolor="none")
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

axinset1.xaxis.set_minor_locator(AutoMinorLocator())
axinset1.yaxis.set_minor_locator(AutoMinorLocator())
axinset2.yaxis.set_minor_locator(AutoMinorLocator())
axinset2.xaxis.set_minor_locator(AutoMinorLocator())

axinset1.set_xlim(T_lower_bound, T_upper_bound+1)
axinset1.set_ylim(0,)
axinset2.set_xlim(T_lower_bound, T_upper_bound+1)
axinset2.set_ylim(0,)


ax2.set_xlim(left=T_lower_bound,right=T_upper_bound+1)
ax2.set_ylim(0,)
ax2.legend(loc="lower left", facecolor="none")
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())

ax3.set_xlim(left=T_lower_bound,right=T_upper_bound+1)
ax3.set_ylim(0,)
ax3.legend(loc="upper left", facecolor="none")
ax3.xaxis.set_minor_locator(AutoMinorLocator())
ax3.yaxis.set_minor_locator(AutoMinorLocator())

ax4.set_xlim(left=T_lower_bound,right=T_upper_bound+1)
ax4.set_ylim(0,)
ax4.legend(loc="upper left", facecolor="none")
ax4.xaxis.set_minor_locator(AutoMinorLocator())
ax4.yaxis.set_minor_locator(AutoMinorLocator())


ax6.axvline(x=873, color="black", ls="dashed")
ax6.axvline(x=1173, color="black", ls="dashed")
ax6.axvline(x=delta_oxygen_inversion_temperature, color="black", ls="dashed")

ax6.set_xlim(T_lower_bound, T_upper_bound+1)
ax6.set_xlabel("T [K]")
ax6.set_ylabel("$\\delta$")
ax6.xaxis.set_minor_locator(AutoMinorLocator())
ax6.yaxis.set_minor_locator(AutoMinorLocator())
ax6.legend()

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
secyax.set_ylabel("% of initial Sr content $\\frac{100 \\cdot x_{eq}}{x_0}$")

secyax3= ax3.secondary_yaxis("right", functions=(yaxconvert, yaxinvert))
secyax3.set_ylabel("% of initial Sr content $\\frac{100 \\cdot x_{eq}}{x_0}$")

secyax4= ax4.secondary_yaxis("right", functions=(yaxconvert, yaxinvert))
secyax4.set_ylabel("% of initial Sr content $\\frac{100 \\cdot x_{eq}}{x_0}$")

secyax.yaxis.set_minor_locator(AutoMinorLocator())
secyax3.yaxis.set_minor_locator(AutoMinorLocator())
secyax4.yaxis.set_minor_locator(AutoMinorLocator())

plt.show()

fig.savefig(local_path + "figs/dry_air_conditions.svg", dpi=300, transparent=True)
fig2.savefig(local_path + "figs/dry_air_condition_deltaG.svg", dpi=300, transparent=True)
#fig3.savefig(local_path + "figs/dry_air_cond_po2_dep.svg", dpi=300, transparent=True)
#fig4.savefig(local_path + "figs/dry_air_cond_p_dep.svg", dpi=300, transparent=True)


