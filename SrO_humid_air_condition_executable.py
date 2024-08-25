import matplotlib.pyplot as plt
from lib.SrO_humid_models import *

x = 0.4 #molar fraction of Sr
x_O2 = 0.21
x_H2O = 0.08
P = 1 #atm
T_lower_bound = 700
T_upper_bound = 1000
T_range = np.arange(T_lower_bound,T_upper_bound,1) #K
#numpy imported chemical potentials, which is imported in humid models

p_O2 = x_O2 * P
p_H2O = x_H2O * P

fig,ax = plt.subplots(layout='constrained')


V_Sr = case1(T_range, x, p_O2, p_H2O, P=1)
#keep in mind that the pH2 is determined inside the case1 function by pO2 and pH2O by considering the equilibrium between these three gases.
ax.plot(T_range, V_Sr, label="H$_{2(g)}$")

V_Sr = case2(np.arange(400,1400), x, p_H2O, P=1)
ax.plot(np.arange(400,1400), V_Sr, label="[2H]$_{La}^{'}$")

V_Sr = case3(T_range, x=0.4, p_O2 = 0.21, p_H2O = 0.08, P=1)
ax.plot(T_range, V_Sr, label="$\\frac{1}{2}$H$_{2(g)}$+[H]$_{La}^{''}$")


V_Sr = case4(T_range, x, p_O2, p_H2O, P=1)
ax.plot(T_range, V_Sr, label="H$_{2(g)}$, Sr$_{bulk}$")

#===================================================
#PLOTTING OPTIONS
ax.set_title("Humid conditions")
ax.set_xlabel("T[K]")
ax.set_ylabel("[V\'\'\'$_{La}$]")


ax.set_xlim(left=400,right=T_upper_bound)
ax.set_ylim(0,)
ax.legend(loc="upper right")

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
