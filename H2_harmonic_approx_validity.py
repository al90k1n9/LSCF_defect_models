import numpy as np
from lib.chemical_potentials import chem_pot_H2, ev2J_p_mol, N_avagadro, R
from lib.dft_energies_0K import E_DFT_H2, E_DFT_H2O
import matplotlib.pyplot as plt
from lib.auxilliary_functions import pH2_giver

ev2J = 1.60219*10**(-19)
kB = 1.380649 * 10**(-23) #J/K
hbar = 1.054571817*10**(-34) #reduced planck's constant in J.s
T_range = np.arange(700,1400)
beta = 1/(kB*T_range)


k = 31.7 #eV/A^2
print("spring constant k in eV/A^2: ", k)
k = k*1e20 * ev2J #J/m^2

mass_hydrogen = 1.67e-27 #kg
reduced_mass = mass_hydrogen**2/(2*mass_hydrogen)
omega = np.sqrt(2*k/reduced_mass)
f = omega/(2*np.pi)
wave_number = f/299792458



print("spring constant k of H2 in J/m^2",k)
print("reduced mass of H2 in kg: ", reduced_mass)
print("harmonic frequncy w0 in m-1: ",wave_number)

#print("order n,\tn-th term in Z,\t\tcumulative sum of Z,\t\tF = -kb*T*ln(Z) in eV")
Z = 0
for n in range(0,4):
    Z += np.exp(-hbar * beta*omega * (n+0.5))
    #F = -kB * T_range * np.log(Z) * N_avagadro # in J/mol
    #print(n,"\t\t",np.exp(-hbar * beta*omega * (n+0.5)),"\t\t",Z,"\t\t", F/ev2J_p_mol)

F = -kB * T_range * np.log(Z) * N_avagadro # in J/mol

delta_mu_H2 = chem_pot_H2(T_range, E_DFT_H2) - E_DFT_H2
#partial_pressure_term = []
#for T in T_range:
#    partial_pressure_term.append(R*T*np.log(pH2_giver(T, 0.08, 0.21)))
#partial_pressure_term = np.asarray(partial_pressure_term)
#
#delta_mu_H2 += partial_pressure_term


print("Hé_energy at OK: ",E_DFT_H2/ev2J_p_mol)
fig, ax = plt.subplots(layout="constrained")
ax.plot(T_range, F/ev2J_p_mol, label="F$_{vib, H2}$")
ax.plot(T_range, delta_mu_H2/ev2J_p_mol, label="$\Delta \mu_{, H2}$")
ax.legend()
ax.set_xlabel("T[K]")
ax.set_ylabel("Energy[eV]")
plt.show()
