
from lib.dft_energies_0K import *
from lib.chemical_potentials import *
from lib.auxilliary_functions import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

#N_avagadro, ev2j, ev2J_p_mol defined in dft_energies_0k
x_H2O = 0.08
x_O2 = 0.21
P = 1 #Bar
kb = 1.380649 * 10**(-23) #J/K
m_H2O = 18.01528 / (N_avagadro*1000) #in kg
m_O2 = 32 / (N_avagadro*1000) #in kg
hbar = 1.054571817*10**(-34) #reduced planck's constant in J.s


oxygen_adsorption = -0.6 #eV
oxygen_adsorption *= ev2J_p_mol #J/mol

print("oxygen adsorption ", oxygen_adsorption*N_avagadro)
print("boltzmann constant in J ", kb)
print("adsorption energy in eV ", E_ads/ev2J_p_mol)


def surface_coverage_modified(T, x_H2O, E_ads, P):
    chemical_potential = cp_H2O(T, E_DFT_H2O=0, P=x_H2O*P) +  zpe_H2O #J/mol
    #BE CAREFUL YOU NEED TO GET THE CHEMICAL POTENTIAL CORRECTION; IT'S IMPORTANT TO PUT THE DFT ENERGY PARAMETER TO ZERO
    exp_term = np.exp(-(-E_ads+chemical_potential)/(R*T))
    theta = 1/(1+exp_term)
    return theta


def langumuir_isotherm(T, x_H2O, E_ads, P):
    chemical_potential = cp_H2O(T, E_DFT_H2O=0, P=x_H2O*P) +  zpe_H2O #J/mol
T_range = np.arange(600, 1301)


theta_list_def_model = surface_coverage_H2O(T_range, x_H2O, E_ads, P=1)
theta_modified = surface_coverage_modified(T_range, x_H2O, E_ads, P=1)




fig,ax = plt.subplots(layout="constrained")


ax.plot(T_range, theta_list_def_model, label="in def model")
ax.plot(T_range, theta_modified, label="modified")


#x_H2O_range = np.linspace(0.02,1,5)
#for x_H2O in x_H2O_range:
#    theta_list = surface_coverage_H2O(T_range, x_H2O, E_ads, P=1)
#    ax.plot(T_range, theta_list_def_model, label="$x_{H_2O} = $" + str(x_H2O))

ax.set_xlim(T_range[0], T_range[-1])
ax.set_ylim(0,1.01)

ax.set_xlabel("T [K]")
ax.set_ylabel("$\\theta_{H_2O}$")

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
#ax.set_xticks([min(T_range), max(T_range)])
#ax.set_xticklabels([f'{min(T_range)}', f'{max(T_range)}'])
ax.legend()

fig.savefig("surf_coverage.svg", format="svg", dpi=300, transparent=True)
plt.show()
