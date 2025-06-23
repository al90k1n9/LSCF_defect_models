"""All results related to adsorption models are presented here.
"""

from lib.dft_energies_0K import *
from lib.chemical_potentials import *
from lib.auxilliary_functions import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
from tqdm import tqdm
import os

local_path = os.path.dirname(os.path.abspath(__file__))
local_path += "/"

default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

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

print("oxygen adsorption ", oxygen_adsorption/ev2J_p_mol)
print("boltzmann constant in J ", kb)
print("adsorption energy in eV ", E_ads/ev2J_p_mol)


def comp_ads(T, x_H2O, x_O2, E_ads_H2O, E_ads_O2, P):
    """Competitive adsorption model

    Parameters
    ----------
    T : float or numpy array
        Temperature or temperature window in Kelvins.
    x_H2O : float
        Water vapour partial pressure
    x_O2 : float
        Oxygen gas partial pressure
    E_ads_H2O : float
        Water molecule energy in J/mol from DFT
    E_ads_O2 : float
        Oxygen molecule energy in J/mol from DFT
    P : float
        Total pressure that defaults to the ambient total pressure of 1 atm/Bar

    Returns
    -------
    tuple of length 4
        Returns a tuple of four elements of same type as T: oxygen surface coverage, water surface coverage, oxygen gas chemical potential in J/mol, water vapour chemical potential in J/mol.

    """
    chemical_potential_O2 = chem_pot_O2(T, E_DFT_O2 = 0, P=x_O2*P)
    chemical_potential_H2O = chem_pot_H2O(T, E_DFT_H2O = 0, P=x_H2O*P)
    exp_term_H2O = np.exp(-(E_ads_H2O - chemical_potential_H2O)/(R*T))
    exp_term_O2 = np.exp(-(E_ads_O2 - chemical_potential_O2)/(R*T))
    theta_O2 = 1/(1+1/exp_term_O2 + exp_term_H2O/exp_term_O2)
    theta_H2O = 1/(1+1/exp_term_H2O + exp_term_O2/exp_term_H2O)
    return (theta_O2, theta_H2O, chemical_potential_O2, chemical_potential_H2O)

T_range = np.arange(400, 1301, 0.1)


inversion_temp_list=[]
half_coverage_temp_list=[]

theta_list_def_model, chemical_potential = surface_coverage_H2O(T_range, x_H2O=0.08, E_ads=E_ads, P=1, chem_pot = 0)
fig0, ax0 = plt.subplots(layout="constrained")
ax0.plot(T_range, theta_list_def_model)
Delta_G = E_ads - chem_pot_H2O(T_range,E_DFT_H2O = 0, P=x_H2O) - zpe_H2O
energy_ax = ax0.twinx()

energy_ax.plot(T_range, Delta_G/ev2J_p_mol, color="black")
energy_ax.axhline(y=0, color="black", linestyle="dashed")


fig,ax = plt.subplots(layout="constrained")
for x_H2O in [0.08, 0.2, 0.5]:
	theta_list_def_model, chemical_potential = surface_coverage_H2O(T_range, x_H2O, E_ads, P=1)
	Delta_G = E_ads - chem_pot_H2O(T_range,E_DFT_H2O = 0, P=x_H2O) - zpe_H2O
	#energy_axes = ax.twinx()
	ax.plot(T_range, theta_list_def_model, label="$x_{H_2O}=$"+str(x_H2O))


x_H2O_range = np.linspace(0.01, 0.99, 200)
for x_H2O in x_H2O_range:
	theta_list_def_model, chemical_potential = surface_coverage_H2O(T_range, x_H2O, E_ads, P=1)
	Delta_G = E_ads - chem_pot_H2O(T_range,E_DFT_H2O = 0, P=x_H2O) - zpe_H2O
	#energy_axes = ax.twinx()
	#ax.plot(T_range, theta_list_def_model, label="$x_{H_2O}=$"+str(x_H2O))
	for index in range(1,len(T_range)):
		if Delta_G[index-1]<0 and Delta_G[index]>=0: inversion_temp_list.append(T_range[index])
		if theta_list_def_model[index-1]>0.5 and theta_list_def_model[index]<=0.5: half_coverage_temp_list.append(T_range[index])

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
ax.legend(facecolor="none")

ax0.set_xlim(T_range[0], T_range[-1])
ax0.set_ylim(0,1.01)

ax0.set_xlabel("T [K]")
ax0.set_ylabel("$\\theta_{H_2O}$")

ax0.xaxis.set_minor_locator(AutoMinorLocator())
ax0.yaxis.set_minor_locator(AutoMinorLocator())

energy_ax.set_ylabel("$\Delta_r^{ads} G(T,p=p_{H_2O})$ [eV]")
energy_ax.yaxis.set_minor_locator(AutoMinorLocator())


#ax.set_xticks([min(T_range), max(T_range)])
#ax.set_xticklabels([f'{min(T_range)}', f'{max(T_range)}'])

fig2, ax2 = plt.subplots(layout="constrained")
ax2.plot(inversion_temp_list, half_coverage_temp_list, marker="s", linestyle="")
ax2.plot(inversion_temp_list, inversion_temp_list, linestyle="dashed", color="black", label="x=y")

ax2.legend(facecolor="none")

ax2.set_xlabel("Inversion temperature [K]")
ax2.set_ylabel("Half coverage temperature [K]")

ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())

fig3, ax3 = plt.subplots(layout="constrained")

ax3.plot(x_H2O_range, inversion_temp_list)

ax3.set_xlabel("$x_{H_2O}$")
ax3.set_ylabel("Inversion temperature [K]")

ax3.set_xlim(0.01, 0.1)
ax3.set_ylim(1075, 1165)

ax3.xaxis.set_minor_locator(AutoMinorLocator())
ax3.yaxis.set_minor_locator(AutoMinorLocator())


ax3inset = ax3.inset_axes([0.45, 0.1, 0.5, 0.5])
ax3inset.set_xlim(0,1)
ax3inset.plot(x_H2O_range, inversion_temp_list)


ax3inset.set_facecolor("none")
ax3inset.xaxis.set_minor_locator(AutoMinorLocator())
ax3inset.yaxis.set_minor_locator(AutoMinorLocator())


fig4, ax4 = plt.subplots(layout="constrained")
theta_O2, theta_H2O, chemical_potential_O2, chemical_potential_H2O = comp_ads(T_range, 0.08, 0.21, E_ads, oxygen_adsorption, 1)

ax4.plot(T_range, theta_H2O)
ax4twinax = ax4.twinx()
ax4twinax.plot(T_range, theta_O2, color=default_colors[1])

#ax4.axvline(x=940, color="black", linestyle="dashed")

ax4.set_xlabel("T [K]")
ax4.set_ylabel("$\\theta_{H_2O}$")


ax4twinax.set_ylabel("$\\theta_{O_2}$")



ax4twinax.yaxis.set_minor_locator(AutoMinorLocator())
ax4.xaxis.set_minor_locator(AutoMinorLocator())
ax4.yaxis.set_minor_locator(AutoMinorLocator())

ax4.set_xlim(T_range[0], T_range[-1])
ax4.set_ylim(0,)
ax4twinax.set_ylim(0,)

fig5, ax5 = plt.subplots(layout="constrained")

ax5.plot(T_range, (E_ads-chemical_potential_H2O)/ev2J_p_mol, label="$\Delta_r^{C5}G^*(T, p=p_{H_2O})$")
ax5.plot(T_range, (oxygen_adsorption - chemical_potential_O2)/ev2J_p_mol, label="$\Delta_r^{C4}G^*(T, p=p_{O_2})$")
difference = ((oxygen_adsorption - chemical_potential_O2)-(E_ads-chemical_potential_H2O))/ev2J_p_mol

#ax5.plot(T_range, np.abs(difference), color="black", label="|$\Delta_r^{C4}G^*(T, p=p_{O_2})-\Delta_r^{C5}G^*(T, p=p_{H_2O})$|")



ax5.set_xlabel("T [K]")
ax5.set_ylabel("[eV]")

ax5.legend(facecolor="none")

ax5.xaxis.set_minor_locator(AutoMinorLocator())
ax5.yaxis.set_minor_locator(AutoMinorLocator())

ax5.set_xlim(T_range[0], T_range[-1])
item = 1

#for oxygen_adsorption in tqdm(np.linspace(-3.0, 0, 1000)):
#	oxygen_adsorption *= ev2J_p_mol
#	theta_O2, theta_H2O, chemical_potential_O2, chemical_potential_H2O = comp_ads(T_range, 0.08, 0.21, E_ads, oxygen_adsorption, 1)
#	fig6, ax6 = plt.subplots(layout="constrained")
#	fig6.get_layout_engine().set(wspace=0.3, hspace=0.3)
#
#	ax6.plot(T_range, theta_H2O)
#	ax6_alternate = ax6.twinx()
#	ax6_alternate.plot(T_range, theta_O2, color=default_colors[1])
#
#	ax6.set_xlabel("T[K]")
#	ax6.set_ylabel("$\\theta_{H_2O}$")
#	ax6_alternate.set_ylabel("$\\theta_{O_2}$")
#	ax6.set_xlim(T_range[0], T_range[-1])
#	ax6.set_ylim(0,)
#	ax6_alternate.set_ylim(0,)
#
#	ax6.xaxis.set_minor_locator(AutoMinorLocator())
#	ax6.yaxis.set_minor_locator(AutoMinorLocator())
#	ax6_alternate.yaxis.set_minor_locator(AutoMinorLocator())
#	ax6.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
#	ax6_alternate.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
#
#	ax6.set_title("E$_{O_2}^{ads}$ = " + str(round(oxygen_adsorption/ev2J_p_mol,2)) + " eV")
#	fig6.savefig(local_path + "figs/gif_comp_surf_coverage/"+str(item)+".png", dpi=300, format="png", transparent=True)
#	item += 1
#	plt.close()
#fig0.savefig(local_path + "figs/surface_coverage.png", dpi=300, transparent=True, format="png")
#fig.savefig(local_path + "figs/surface_coverage_xH2O.png", dpi=300, transparent=True, format="png")
#fig2.savefig(local_path + "figs/inversion_temp_half_coverage_temp.png", dpi=300, transparent=True, format="png")
fig3.savefig(local_path + "figs/inverstion_temp_xH2O.svg", dpi=300, transparent=True, format="png")
#fig4.savefig(local_path + "figs/comp_adsorption.png", format="png", dpi=300, transparent=True)
#fig5.savefig(local_path + "figs/comp_adsorption_delta_Gs.png", format="png", dpi=300, transparent=True)
plt.show()
