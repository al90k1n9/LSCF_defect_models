import matplotlib.pyplot as plt 
from lib.SrO_dry_air_models import *
#chemical potentials are imported in models
#numpy imported chemical potentials
from matplotlib.ticker import  AutoMinorLocator

x = 0.4 #molar fraction of Sr
x_O2 = 0.21
P = 1 #atm
T_lower_bound = 700
T_upper_bound = 1400
T_range = np.arange(T_lower_bound,T_upper_bound,1) #K

p_O2 = x_O2 * P


fig,ax = plt.subplots(layout='constrained')
fig2, ax2 = plt.subplots(layout ="constrained")

ylist, delta_G_list = case1(T_range, x, p_O2, P)
plt_element_case1 = ax.plot(T_range, ylist, label ="case 1")
plt_element_case1_dG = ax2.plot(T_range, np.asarray(delta_G_list)/ev2J_p_mol, label="case 1")

ylist, delta_G_list = case2(T_range, x, p_O2, P)
plt_element_case2, = ax.plot(T_range, ylist, label = "case 2")
plt_element_case2_dG = ax2.plot(T_range, np.asarray(delta_G_list)/ev2J_p_mol, label="case 2")

ylist, delta_G_list = case3(T_range, x=0.4)
plt_element_case3, =ax.plot(T_range, ylist, label ="case 3")
plt_element_case3_dG = ax2.plot(T_range, np.asarray(delta_G_list)/ev2J_p_mol, label="case 3")

ylist, delta_G_list = case4(T_range, x = 0.4)
plt_element_case4, = ax.plot(T_range, ylist, label ="case 4")
plt_element_case4_dG = ax2.plot(T_range, np.asarray(delta_G_list)/ev2J_p_mol, label="case 4")

ylist, delta_G_list = case5(T_range, x, p_O2, P)
plt_element_case5 = ax.plot(T_range, ylist, label ="case 5")
plt_element_case5_dG = ax2.plot(T_range, np.asarray(delta_G_list)/ev2J_p_mol, label="case 5")

ylist, delta_G_list = case6(T_range, x=0.4)
plt_element_case6 = ax.plot(T_range, ylist, label ="case 6")
plt_element_case6_dG = ax2.plot(T_range, np.asarray(delta_G_list)/ev2J_p_mol, label="case 6")

ax.set_title("Dry air conditions")
ax.set_xlabel("T[K]")
ax.set_ylabel("[V\'\'\'$_{La}$]")


#ax2.set_title("Dry air conditions")
ax2.set_xlabel("T[K]")
ax2.set_ylabel("${\Delta}G$[eV]")


ax.set_xlim(left=T_lower_bound,right=T_upper_bound)
ax.set_ylim(0,)
ax.legend(loc="upper left")
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

ax2.set_xlim(left=T_lower_bound,right=T_upper_bound)
ax2.set_ylim(0,)
ax2.legend(loc="lower right")
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())


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



