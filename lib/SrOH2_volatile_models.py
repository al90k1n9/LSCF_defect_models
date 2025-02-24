import sys
sys.path.append("H:\Documents\SrO_defect_model")


from lib.chemical_potentials import *
from lib.dft_energies_0K import * 
from lib.auxilliary_functions import * 
import matplotlib.pyplot as plt


data = np.genfromtxt("./lib/sroh2_factsage.csv", delimiter=";")

T_range= data[:,0]

mu_H2O = cp_H2O(T_range, E_DFT_H2O = E_DFT_H2O)
mu_SrOH2 = data[:,2] + mu_H2O + E_SrO - 0.5*ev2J_p_mol

delta_G = (2* mu_SrOH2 + E_LSCF_slab_Sr_surf_O_sub_surf-2*mu_H2O-E_LSCF_slab)/2
V_Sr = np.zeros(len(T_range))

for index in range(len(T_range)):
    A = 1
    B = np.exp(-delta_G[index]/(R*T_range[index]))
    C = -0.4*np.exp(-delta_G[index]/(R*T_range[index]))
    V_Sr[index] = quadratic_model(A,B,C)

fig, ax = plt.subplots(layout="constrained")

ax.plot(T_range, delta_G/ev2J_p_mol)
ax.set_xlabel("T [K]")
ax.set_ylabel("delta_G [eV]")


fig2, ax2 = plt.subplots(layout="constrained")
ax2.plot(T_range, V_Sr)

ax2.set_xlabel("T[K]")
ax2.set_ylabel("V_Sr")

plt.show()