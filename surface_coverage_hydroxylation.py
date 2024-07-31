
from lib.dft_energies_0K import *
import numpy as np
import matplotlib.pyplot as plt


#N_avagadro, ev2j, ev2J_p_mol defined in dft_energies_0k
kB = 1.380649 * 10**(-23) #J/K
m_H2O = 18.01528 / (N_avagadro*1000) #in kg
hbar = 1.054571817*10**(-34) #reduced planck's constant in J.s
x_H2O = 0.08
P = 1 #Bar

E_ads = (E_LSCF_hydroxilated - E_LSCF_slab - 2*E_DFT_H2O)/2 * ev2J/ev2J_p_mol #J
print(E_ads/ev2J) #eV

T_range = np.arange(400, 2000)
theta = []
for T in T_range:
    p_zero = kB*T*(m_H2O * kB * T/(2*np.pi * hbar**2))**(3/2) * np.exp(E_ads/(kB*T))
    coverage = x_H2O*P*1e5/(x_H2O*P*1e5 + p_zero)
    theta.append(coverage)


plt.plot(T_range, theta)
plt.ylabel("$\\theta$")
plt.xlabel("T[K]")
plt.ylim(-0.1,1.1)
plt.xlim(400,2000)
plt.title("x$_{H_2O}$ = " + str(x_H2O))

plt.show()
