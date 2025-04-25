from lib.dft_energies_0K import *
from lib.chemical_potentials import R
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FuncFormatter

def case2(T, theta, x=0.4, p_H2O = 0.08, P=1):
    delta_E = (E_LSCF_double_hydrogenated + 2*E_SrO_epitax - E_LSCF_hydroxilated )/2 + E_int
    #print(delta_E/ev2J_p_mol, " of case 2")
    delta_G = delta_E
    N = theta/(1-theta) * np.exp(-delta_G/(R*T))
    return N/(1+N)*x


T = 1000
theta_range = [0]
V_Sr=[0]
threshold_decomp = 0.3

while V_Sr[-1]<threshold_decomp:
    theta_step = (1-theta_range[-1])/4
    new_theta = theta_range[-1] + theta_step
    theta_range.append(new_theta)
    V_Sr.append(case2(T, new_theta))


fig,ax = plt.subplots(layout="constrained")
#axinset = ax.inset_axes([0.1, 0.2, 0.5,0.5])

ax.plot(1-np.asarray(theta_range), V_Sr)
#axinset.plot(theta_range, V_Sr)

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

ax.set_xlabel("$\\theta_{H_2O}$")
ax.set_ylabel("$x_{eq} = [V_{Sr}\'\'\']$")


ax.set_ylim(0,)
#ax.set_xlim(1, 0)

ax.set_xscale("log")


def convert_ticks(x, pos):
    # For example, multiply the tick value by 10 and format as integer
    return f'{(-x+1)}'

# Apply the function to format the x-axis ticks
ax.xaxis.set_major_formatter(FuncFormatter(convert_ticks))

ax.invert_xaxis()

#axinset.xaxis.set_minor_locator(AutoMinorLocator())
#axinset.yaxis.set_minor_locator(AutoMinorLocator())

#axinset.set_ylim(0,)
#axinset.set_xlim(lower_bound,1)

fig.savefig("figs/surf_cov_sensitivity.svg", format="svg", dpi=300, transparent=True)


plt.show()

