from lib.dft_energies_0K import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import FormatStrFormatter as fsf
from matplotlib.ticker import AutoMinorLocator

acell = 1.4641598078E+01 
Bohr2Angstrom = 0.5291777
Ha2eV = 27.2114
ev2J = 1.60219*10**(-19)
kB = 1.380649 * 10**(-23) #J/K
hbar = 1.054571817*10**(-34) #reduced planck's constant in J.s
T_range = np.arange(0,1400)

O58 = np.asarray([7.4076004134E-01,  2.1591012494E-01,  4.1481175529E-01])
O87 = np.asarray([8.1539519077E-01,  3.1418832058E-01,  4.1970138000E-01])
O82 = np.asarray([7.5741977865E-01,  1.9588505177E-01,  9.8492142701E-01])
O88 = np.asarray([8.2533191524E-01,  2.9543639395E-01,  9.7782933190E-01])

first_bond = O87 - O58
first_bond_length = np.sqrt(np.matmul(np.transpose(first_bond), first_bond)) * acell * Bohr2Angstrom

second_bond = O88 - O82
second_bond_length = np.sqrt(np.matmul(np.transpose(second_bond), second_bond)) *acell * Bohr2Angstrom

bond_length = (first_bond_length + second_bond_length)/2

delta_bond_lengths_percent = [0, 15, -15, -10, -5, 10, 5]
bond_lengths_angstrom = []
delta_bond_lengths_angstrom = np.asarray(delta_bond_lengths_percent)/100 * bond_length
bond_lengths_angstrom = delta_bond_lengths_angstrom + bond_length 
single_hydrogenation_configs = single_hydrogenation_configs/ev2J_p_mol

#shifting energies
single_hydrogenation_configs -= single_hydrogenation_configs[0]

def quadratic_fit(x, k):
    return 1/2 * k * x**2


force_constant_k, cov = curve_fit(quadratic_fit, delta_bond_lengths_angstrom, single_hydrogenation_configs, [0])

#k = 31.7 #eV/A^2
k = force_constant_k*1e20 * ev2J #J/m^2


mass_hydrogen = 1.67e-27 #kg
mass_oxygen = 2.6566962 * 10**-26 #kg
reduced_mass = (mass_hydrogen*mass_oxygen)/(mass_oxygen + mass_hydrogen)
#reduced_mass = mass_hydrogen

omega = np.sqrt(k/reduced_mass)
f = omega/(2*np.pi) #frequency charactersitic
wave_number = f/299792458




def vibrational_correction_term(T):
    Z= 0
    beta = 1/(kB*T)
    #print("order n\t\t Z \t\t\t nth term \t\t F")
    print(T)
    for n in range(0,4):
        Z += np.exp(-hbar * beta*omega * (n+0.5))
        #F_vib = -kB * T * np.log(Z)
        #print(n,"\t\t", Z,"\t\t", np.exp(-hbar * beta*omega * (n+0.5)),"\t", F_vib/ev2J)

    F_vib = -kB * T * np.log(Z)*N_avagadro #for the sake of consistency keeping every quantity in Joules per mol units
    return F_vib




if __name__ == "__main__":
    print("first OH bond length in A", first_bond_length)
    print("second OH bond length in A", second_bond_length)
    print("average of the above two", bond_length, "in A")

    print("force constant k", force_constant_k, "eV/A^2")
    fig1, ax1 = plt.subplots(layout="constrained")
    ax1.plot(delta_bond_lengths_angstrom, single_hydrogenation_configs, marker="s", ls="")

    fit_xlist = np.linspace(min(delta_bond_lengths_angstrom),max(delta_bond_lengths_angstrom), 1000)
    

    ax1.plot(fit_xlist, quadratic_fit(fit_xlist, force_constant_k[0]), color = "black", ls="dashed")
    ax1.set_xlabel("u [A]")
    ax1.set_ylabel("delta total energy [eV]")
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.set_ylim(0,)

    print("\n\nDetermination of vibrational correction term to the chemical potential")
    print("spring constant k of OH in eV/A^2", force_constant_k)
    print("spring constant k of OH in J/m^2",k)
    print("reduced mass of OH in kg: ", reduced_mass)
    print("harmonic frequncy w0 in m-1: ",wave_number)
    #print("F correction to the chemical potential in eV: ", F_vib/ev2J_p_mol)
    print("verifiying the function")
    print("vibrational correction term at T =1000K in (J/mol):",vibrational_correction_term(1000))
    F_vib = vibrational_correction_term(T_range)
    fig2, ax2 = plt.subplots(layout="constrained")
    ax2.plot(T_range, F_vib/ev2J_p_mol)
    ax2.set_xlabel("T[K]")
    ax2.set_ylabel("F$_{O-H}^{vib}$ [eV]")
    ax2.set_xlim(min(T_range),max(T_range)+1)
    #ax2.yaxis.set_major_formatter(fsf('%.2e'))
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())


    fig1.savefig("figs/H_vibrations_parabolic_fit.svg", format="svg", transparent=True, dpi=300)
    fig2.savefig("figs/H_vibrations_F_vib.svg", format="svg", transparent=True, dpi=300)
    plt.show()
