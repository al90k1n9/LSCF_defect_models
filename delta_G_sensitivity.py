import matplotlib.pyplot as plt
from lib.SrO_hydroxylated_models import case2 as hydrox
from lib.SrO_humid_models import case2 as humid
from lib.SrO_humid_models import case1 as humid_standard
from lib.SrO_dry_air_models import *
from matplotlib.ticker import  AutoMinorLocator
import os

local_path = os.path.dirname(os.path.abspath(__file__))
local_path += "/"

x_O2 = 0.21
P = 1 #atm
T_lower_bound = 600 
T_upper_bound = 1200
T_range = np.arange(T_lower_bound,T_upper_bound,1) #K
x = 0.4

fig, ax = plt.subplots(layout="constrained")
fig2, ax2 = plt.subplots(layout="constrained")
fig3, ax3 = plt.subplots(layout="constrained")
fig4, ax4 = plt.subplots(layout="constrained")
#ax4inset = ax4.inset_axes([0.45,0.175,0.5,0.5], facecolor="none")
fig5, ax5 = plt.subplots(layout="constrained")


delta_oxygen_parameters = delta_oxygen_interpolater(plot=0)
#delta_oxygen_parameters = np.asarray([0,0])
delta_oxygen_inversion_temperature = -delta_oxygen_parameters[1]/delta_oxygen_parameters[0]

ylist, delta_G_list, delta_oxygen_list = case4(T_range, x = 0.4, delta_oxygen_parameters=delta_oxygen_parameters)
ax.plot(T_range, ylist, label ="shift = 0", color="black", ls="dashed")
ax2.plot(T_range, np.asarray(delta_G_list)/ev2J_p_mol, label="shift = 0", color="black", ls="dashed")

for shift in -np.arange(0.1, 1, 0.2):
    shift_value = shift * ev2J_p_mol
    ylist, delta_G_list, delta_oxygen_list = case4(T_range, x = 0.4, delta_oxygen_parameters=delta_oxygen_parameters, sensitivity_shift=shift_value)
    ax.plot(T_range, ylist, label ="shift = " + str(round(shift, 2)) + " ev")
    ax2.plot(T_range, np.asarray(delta_G_list)/ev2J_p_mol, label="shift = " + str(round(shift, 2)) + " ev")

V_Sr, delta_G, theta_list =  hydrox(T_range)
ax3.plot(T_range, V_Sr, color="black", ls="dashed", label="shift=0")


for shift in -np.arange(0.1, 1, 0.2):
    shift_value = shift * ev2J_p_mol
    V_Sr, delta_G, theta_list =  hydrox(T_range, sensitivity_shift = shift_value)
    ax3.plot(T_range, V_Sr, label="shift = " + str(round(shift, 2)) + " ev")

fix_T = [1000]
ylist, delta_G_list, delta_oxygen_list = case4([1000], x = 0.4)
dry_case = [ylist[0]]
V_Sr, delta_G, theta_list =  hydrox([1000])
hydrox_case = [V_Sr[0]]
V_Sr, delta_G_range = humid([1000])
humid_case = [V_Sr[0]]
V_Sr, delta_G_range, ph2_list= humid_standard([1000])
humid_standard_case = [V_Sr[0]]

shift_step = 0.01
shift_range = -np.arange(shift_step, 1, shift_step)
for shift in shift_range:
    shift_value = shift *ev2J_p_mol

    V_Sr, delta_G, delta_oxygen_list = case4(fix_T, sensitivity_shift= shift_value)
    dry_case.append(V_Sr[0])

    V_Sr, delta_G, theta_list = hydrox(fix_T, sensitivity_shift= shift_value)
    hydrox_case.append(V_Sr[0])

    V_Sr, delta_G_range = humid(fix_T, sensitivity_shift = shift_value)
    humid_case.append(V_Sr[0])

    V_Sr, delta_G_range, ph2_list= humid_standard(fix_T, sensitivity_shift = shift_value)
    humid_standard_case.append(V_Sr[0])

shift_range = np.insert(shift_range, 0, 0)

ax4.plot(shift_range, dry_case, label="Dry-air C4")
#ax4.plot(shift_range, humid_standard_case, label="Humid R3.1")
ax4.plot(shift_range, hydrox_case, label="Humid R3.2")
ax4.plot(shift_range, humid_case, label="Humid R3.3")


#ax4inset.plot(shift_range, dry_case, label="Dry-air C4")
#ax4inset.plot(shift_range, humid_standard_case, label="Humid R3.1")
#ax4inset.plot(shift_range, hydrox_case, label="Humid R3.2")
#ax4inset.plot(shift_range, humid_case, label="Humid R3.3")




Sp_LSCF_CGO_initial = 6659953.51689531 #1/m
a_LSCF = acell_LSCF_slab/2 #m
print(a_LSCF, " a_LSCF")
print(N_avagadro, " avagadro")
Sp_decrease_percent =np.arange(0,100, 0.01)
slope = 2* Sp_LSCF_CGO_initial *a_LSCF/ (100) /0.4 * 1000
print(slope, " slope")
required_xeq = Sp_decrease_percent * slope 

ax5.plot(Sp_decrease_percent, required_xeq)
#ax5.plot(Sp_decrease_percent, 2*required_xeq)

ax5.set_xlim(0,100)
ax5.set_ylim(0,)
ax5.set_xlabel("Decrease in $S_p^{LSCF/pores}$ in percentages")
ax5.set_ylabel("Required Sr release in percentages")
ax5.yaxis.set_minor_locator(AutoMinorLocator())
ax5.xaxis.set_minor_locator(AutoMinorLocator())

ax.set_xlabel("T")
ax.set_ylabel("$[V\'\'\'_{La}]_{eq}=x_{eq}$")

ax2.set_xlabel("T[K]")
ax2.set_ylabel("${\Delta}G^*(T,p)$ [eV]")

def yaxconvert(x):
    return x * 100/0.4

def yaxinvert(x):
    return x *0.4/100

secyax = ax.secondary_yaxis("right", functions=(yaxconvert, yaxinvert))
secyax.set_ylabel("% of initial Sr content $\\frac{100 \\cdot x_{eq}}{x_0}$")

ax.set_xlim(left=T_lower_bound,right=T_upper_bound+1)
ax.set_ylim(0,)
ax.legend(loc="upper left", facecolor="none")
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

ax2.set_xlim(left=T_lower_bound,right=T_upper_bound+1)
ax2.set_ylim(0,)
ax2.legend(loc="lower left", facecolor="none")
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())

ax3.legend(facecolor="none")
ax3.set_xlim(T_lower_bound, T_upper_bound+1)
ax3.set_ylim(0,)
ax3.set_xlabel("T[K]")
ax3.set_ylabel("$[V\'\'\'_{La}]_{eq}=x_{eq}$")
ax3.xaxis.set_minor_locator(AutoMinorLocator())
ax3.yaxis.set_minor_locator(AutoMinorLocator())

secyax3 = ax3.secondary_yaxis("right", functions=(yaxconvert, yaxinvert))
secyax3.set_ylabel("% of initial Sr content $\\frac{100 \\cdot x_{eq}}{x_0}$")


ax4.legend(facecolor="none")
ax4.set_xlim(min(shift_range), max(shift_range))
ax4.set_ylim(0,)
ax4.set_xlabel("shift energy [eV]")
ax4.set_ylabel("$[V\'\'\'_{La}]_{eq}=x_{eq}$")
ax4.xaxis.set_minor_locator(AutoMinorLocator())
ax4.yaxis.set_minor_locator(AutoMinorLocator())

secyax4 = ax4.secondary_yaxis("right", functions=(yaxconvert, yaxinvert))
secyax4.set_ylabel("% of initial Sr content $\\frac{100 \\cdot x_{eq}}{x_0}$")
secyax4.yaxis.set_minor_locator(AutoMinorLocator())


fig4.savefig(local_path + "figs/delta_G_sensitivity.svg", format = "svg", transparent = True, dpi = 300)
fig5.savefig(local_path + "figs/required_SrO_for_Sp_reduction.svg", format = "svg", transparent = True, dpi = 300)

plt.show()

