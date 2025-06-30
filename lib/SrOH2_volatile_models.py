"""Cases for the formation of SrOH2 volatile
"""

import sys

from lib.chemical_potentials import *
from lib.dft_energies_0K import * 
from lib.auxilliary_functions import * 
import os

local_path = os.path.dirname(os.path.abspath(__file__))
local_path += "/"



def case1(T_range, x=0.4, x_H2O = 0.08, P = 1, volume_fraction_gas= 0.5, delta_oxygen_parameters = [0,0]):
    """Case where atmosphering water vapour is used for the formation of SrOH2 volatile species. LSCF interacts with the gas in the pores.

    Parameters
    ----------
    T_range : numpy array
        Temperature window in Kelvins
    x : float, optional
        Initial Sr content in LSCF, defaults to be standard LSCF composition of 0.4
    x_H2O : float, optional
        Water vapour molar fraction in the gas in pores, defaults to be 0.08 the experimental conditions in Sassone's article
    P : float, optional
        Total pressure, defaults to be ambient gas total pressure of 1 atm
    volume_fraction_gas : float, optional
        Volume fraction of pores in the material. Defaults to be 0.5
    delta_oxygen_parameters : numpy array, optional
        Parameters for linear fitting the oxygen understoichometry as a function of temperature. Array contains slope and offset

    Returns
    -------
    tuple of length 4
        Tuple containing 4 numpy arrays: # Sr vacancies per unit LSCF cell created due to SrOH2 formation, delta_G in J/mol, oxygen understoichimetry in numbers/unit LSCF cell, number of H2O molecules present initially in the pores; all function of temperature.
    
    Warnings
    --------
    Ni_H2O::number of H2O moelcules in pores is a function of temperature as we use ideal gas law.

    """
    p_H2O = x_H2O * P
    V_Sr = np.zeros(len(T_range))
    V_Sr1 = np.zeros(len(T_range))
    V_Sr2 = np.zeros(len(T_range))
    delta_G_list = np.zeros(len(T_range))
    delta_oxygen_list = []
    Ni_H2O_list = []
    for index in range(len(T_range)):
        T = T_range[index]
        mu_H2O = chem_pot_H2O(T, E_DFT_H2O=E_DFT_H2O)
        mu_SrOH2 = chem_pot_SrOH2(T)
        delta_G = float(0.5*(2*mu_SrOH2 + E_LSCF_slab_Sr_surf_O_sub_surf - 2*mu_H2O - E_LSCF_slab))
        delta_G_list[index] = delta_G

        delta_oxygen = delta_oxygen_parameters[0] * T + delta_oxygen_parameters[1]
        if delta_oxygen<0: delta_oxygen = 0 #understoichiometry cannot be negative.
        delta_oxygen_list.append(delta_oxygen)

        Ni_H2O = volume_fraction_gas/(1-volume_fraction_gas) * (acell_LSCF_slab/2)**3/(kb*T) * p_H2O * 1e5
        Ni_H2O_list.append(Ni_H2O)

        a = -1-np.exp(delta_G/(R*T))
        b = Ni_H2O + x + 3
        c = -(x*Ni_H2O + (3-delta_oxygen)*Ni_H2O + (3-delta_oxygen)*x)
        d = (3-delta_oxygen)*x*Ni_H2O
        solutions = cubic_model(a,b,c,d)
        #equation_func = lambda x: equation(x, T=T, p_H2O=  p_H2O, delta_G=delta_G)
        #solution = bisection(0,0.4)

        V_Sr[index] = solutions[0]
    return(np.asarray(V_Sr), np.asarray(delta_G_list), np.asarray(delta_oxygen_list), np.asarray(Ni_H2O_list))


def Ni_H2O_sensitivity(T=1000, x=0.4, volume_fraction_gas= 0.5, delta_oxygen_parameters = [0,0]):
    """Determines the effect of number of water molecules present in the pores initially on the formation of SrOH2 volatile.

    Parameters
    ----------
    T : float
        Temperature in Kelvin
    x : float, optional
        Initial Sr content in LSCF, defaults to standard composition of 0.4
    volume_fraction_gas : float, optional
        Volume fraction of pores in LSCF, defaults to 0.5
    delta_oxygen_parameters : numpy arry, optional
        Linear interpolation parameteres: slopea and offset, to determine oxygen understoichiometry as a function of temperature

    Returns
    -------
    tuple of length 3
        Tuple containing 3 numpy arrays : #H2O molecules initiallyp resent in pores/LSCF cell, #Sr vacancies/LSCF unit cell, molar fraction of H2O in pores; all functions of temperature

    """
    mu_H2O = chem_pot_H2O(T, E_DFT_H2O=E_DFT_H2O)
    mu_SrOH2 = chem_pot_SrOH2(T)
    delta_G = float(0.5*(2*mu_SrOH2 + E_LSCF_slab_Sr_surf_O_sub_surf - 2*mu_H2O - E_LSCF_slab))
    #reference_Ni_H2O = volume_fraction_gas/(1-volume_fraction_gas) * (acell_LSCF_slab/2)**3/(kb*T) * p_H2O * 1e5

    Ni_H2O_list = np.logspace(-7, 0, num=30)
    x_H2O_list = Ni_H2O_list/volume_fraction_gas/((1-volume_fraction_gas) * (acell_LSCF_slab/2)**3/(kb*T)*1e10)
    print(Ni_H2O_list)
    solution_list = []
    for Ni_H2O in Ni_H2O_list:
        a = -1-np.exp(delta_G/(R*T))
        b = Ni_H2O + x + 3
        c = -(x*Ni_H2O + 3*Ni_H2O + 3*x)
        d = 3*x*Ni_H2O
        solution_list.append(cubic_model(a,b,c,d)[0])
    return (np.asarray(Ni_H2O_list), np.asarray(solution_list), np.asarray(x_H2O_list))


if __name__ == "__main__":
    data = np.genfromtxt("./lib/sroh2_factsage.csv", delimiter=";")
    T_range = data[100::100,0]
    T = T_range[99]
    mu_H2O = chem_pot_H2O(T, E_DFT_H2O=E_DFT_H2O)
    mu_SrOH2 = chem_pot_SrOH2(T)
    delta_G = 0.5*(2*mu_SrOH2[0] + E_LSCF_slab_Sr_surf_O_sub_surf - 2*mu_H2O - E_LSCF_slab)

    print("exponential term: ", np.exp(-delta_G/(R*T)))
    print("delta G", delta_G, delta_G/ev2J_p_mol)

    print("temperature T[K]: ",T)

    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator

    fig, ax = plt.subplots(layout="constrained")


    def thermo_constant(x, x_initial = 0.4):
        dinominator = (3-x) * (x_initial-x) * (Ni_H2O- x)
        if dinominator == 0: print(Ni_H2O, x)
        return x**3/dinominator

    volume_fraction_gas = 0.5
    p_H2O = 0.08
    Ni_H2O = volume_fraction_gas/(1-volume_fraction_gas) * (acell_LSCF_slab/2)**3/(kb*T) * p_H2O * 1e5 # #molecules in gas volume to react with one unit LSCF
    print("actual N_H2O = ", Ni_H2O)
    xlist = np.arange(0, Ni_H2O, Ni_H2O/100)
    ylist = []
    polynomial = []
    x_initial = 0.4
    a = -1-np.exp(delta_G/(R*T))
    b = Ni_H2O + x_initial + 3
    c = -(x_initial*Ni_H2O + 3*Ni_H2O + 3*x_initial)
    d = 3*x_initial*Ni_H2O
    for x in xlist:
        ylist.append(thermo_constant(x)-np.exp(-delta_G/(R*T)))

        polynomial.append(a*x**3 + b*x**2 + c*x + d)

    root_function = lambda x:thermo_constant(x,T) - np.exp(-delta_G/(R*T))
    solution = bisection(0, Ni_H2O-1e-12, root_function)
    print("solution from actual function: ", solution)
    print("analytical solution", cubic_model(a,b,c,d))



    ax.plot(xlist, ylist, label = "actual function")
    ax.plot(xlist, polynomial, label = "polynomial")
    ax.axhline(y=0, color="black")
    ax.legend()
    ax.set_xlabel("x")
    ax.set_ylabel("")   
    #ax.set_xlim(0,1e-5)
    ax.set_ylim(-1e-11,1e-11)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    actual_Ni_H2O = Ni_H2O

    Ni_H2O_list = np.logspace(-7, -1, num=30)
    print(Ni_H2O_list)
    solution_list = []
    for Ni_H2O in Ni_H2O_list:
        a = -1-np.exp(delta_G/(R*T))
        b = Ni_H2O + x_initial + 3
        c = -(x_initial*Ni_H2O + 3*Ni_H2O + 3*x_initial)
        d = 3*x_initial*Ni_H2O
        solution_list.append(cubic_model(a,b,c,d)[0])

    fig2, ax2 = plt.subplots(layout="constrained")
    ax2.plot(Ni_H2O_list, solution_list, marker="s")
    ax2.set_xlabel("Ni_H2O [#molecules/gas volume]")
    ax2.set_ylabel("[V$_{Sr}$'''] [#molecules]")
    ax2.axvline(x=actual_Ni_H2O, color="black", label="Ni_H2O = "+str("{:.3e}".format(actual_Ni_H2O)))

    ax2.legend()

    #ax2.set_yscale("log")
    #ax2.set_xscale("log")

    fig2.savefig(local_path + "figs/Ni_H2O_sensitivity.png", dpi=300, transparent=True, format="png")
    plt.show()
