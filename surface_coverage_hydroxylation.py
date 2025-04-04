
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

T_range = np.arange(700, 2000)
theta_list= []
theta_alternative = [] #this list is calculated with the experimental chemical potential from nist-janaf tables
theta_competitive = []

cp_list = []
O2_cp_list = []
simplified_cp_list = []
simplified_O2_cp_list = []
test = []

original_term = []
additional_term = []

oxygen_coverage = []
for T in T_range:
    P0 = kb*T*(m_H2O*kb*T/(2*np.pi*hbar**2))**(3./2)*np.exp(E_ads/(N_avagadro * kb*T))
    simplified_cp = np.log(x_H2O *P* 101325/(kb*T*(m_H2O*kb*T/(2*np.pi*hbar**2))**(3./2)))*kb*T #J
    simplified_O2_cp = np.log(x_O2 *P* 101325/(kb*T*(m_O2*kb*T/(2*np.pi*hbar**2))**(3./2)))*kb*T #J
    theta = x_H2O * P * 101325/(x_H2O * P * 101325 + P0)
    theta_list.append(theta)
    simplified_cp_list.append(simplified_cp)
    simplified_O2_cp_list.append(simplified_O2_cp)


    exp_cp_H2O = cp_H2O(T, E_DFT_H2O=0, P=1) + N_avagadro * kb*T* np.log(x_H2O) #important to have the zero for 
    exp_term = np.exp(-(-E_ads+exp_cp_H2O)/(kb*T*N_avagadro))
    theta_alternative.append(1/(1+exp_term))
    cp_list.append(exp_cp_H2O/N_avagadro)

    exp_cp_O2 = (cp_O2(T, E_DFT_O2=0, P=1) + N_avagadro * kb*T*np.log(x_O2)) #J/mol
    O2_cp_list.append(exp_cp_O2)
    numerator = np.exp((-E_ads + simplified_cp * N_avagadro)/(N_avagadro*kb*T)) * 0.08
    dinominator = 1 + numerator +  0.21 * np.exp((-oxygen_adsorption + simplified_O2_cp * N_avagadro)/(kb*T*N_avagadro))
    theta_competitive.append(numerator/dinominator)
    additional_term.append(np.exp((simplified_O2_cp * N_avagadro)/(kb*T*N_avagadro)))
    original_term.append(numerator)
    oxygen_numerator = 0.21 * np.exp((-oxygen_adsorption + simplified_O2_cp * N_avagadro)/(kb*T*N_avagadro))
    oxygen_coverage.append(oxygen_numerator/dinominator)

theta_list= np.asarray(theta_list)
theta_alternative = np.asarray(theta_alternative)
theta_competitive = np.asarray(theta_competitive)

cp_list = np.asarray(cp_list)
O2_cp_list = np.asarray(O2_cp_list)
simplified_cp_list = np.asarray(simplified_cp_list)
simplified_O2_cp_list = np.asarray(simplified_O2_cp_list)
test = np.asarray(test)

original_term = np.asarray(original_term)
additional_term = np.asarray(additional_term)
oxygen_coverage = np.asarray(oxygen_coverage)



fig,ax = plt.subplots(layout="constrained")
fig2, ax2 = plt.subplots(layout="constrained")
fig3, ax3 = plt.subplots(layout="constrained")
fig4, ax4 = plt.subplots(layout="constrained")


ax.plot(T_range, theta_list, label="simplified chemical potential")
#ax.plot(T_range, theta_alternative, label="exp chemical potential")
#ax.plot(T_range, theta_competitive, label="competetive adsorption")
#ax.plot(T_range, oxygen_coverage, label="oxygen coverage")

ax.set_ylabel("$\\theta$")
ax.set_xlabel("T[K]")
#ax.set_ylim(-0.1,1.1)
ax.set_xlim(T_range[0], T_range[-1])
ax.set_title("Surface coverages \n x$_{H_2O}$ = " + str(x_H2O))

ax.legend(loc="upper right", facecolor="none")

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

ax2.plot(T_range, simplified_cp_list/ev2J, label="translation H2O chemical potential")
ax2.plot(T_range, (simplified_O2_cp_list)/ev2J, label="translation O2 chemical potential")
ax2.plot(T_range, (cp_list)/ev2J, label="exp H2O chemical potential")
ax2.plot(T_range, (O2_cp_list)/ev2J_p_mol, label="exp O2 chemical potential")
ax2.legend(loc="lower left", facecolor="none")
ax2.set_ylabel("chemical potentials [eV]")
ax2.set_xlabel("T[K]")
#ax2.set_ylim(-0.1,1.1)
ax2.set_xlim(T_range[0],T_range[-1])
ax2.set_title("chemical potentials")

ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())

ax3.plot(T_range, additional_term/(original_term))
#ax3.plot(T_range, original_term)

ax3.set_xlabel("T[K]")
ax3.set_ylabel("ration of exp terms")


ax4.plot(T_range, original_term, label="orignal term")
ax4.plot(T_range, additional_term, label="additional term")

ax4.set_xlabel("T[K]")
ax4.set_ylabel("exp terms")
ax4.legend()


fig.savefig("surface_coverage.png", dpi=300, transparent=True)

#simplified cp slope
print("slope of cp in J/K: ", (simplified_cp_list[1]-simplified_cp_list[0])/(T_range[1]-T_range[2]))

plt.show()
