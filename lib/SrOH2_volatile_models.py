import sys

from lib.chemical_potentials import *
from lib.dft_energies_0K import * 
from lib.auxilliary_functions import * 


def case1(T_range, x=0.4, p_H2O = 0.08, P = 1, volume_fraction_gas= 0.5):
    V_Sr = np.zeros(len(T_range))
    V_Sr1 = np.zeros(len(T_range))
    V_Sr2 = np.zeros(len(T_range))
    delta_G_list = np.zeros(len(T_range))

    for index in range(len(T_range)):
        T = T_range[index]
        mu_H2O = cp_H2O(T, E_DFT_H2O=E_DFT_H2O)
        mu_SrOH2 = cp_SrOH2(T)
        delta_G = float(0.5*(2*mu_SrOH2 + E_LSCF_slab_Sr_surf_O_sub_surf - 2*mu_H2O - E_LSCF_slab))
        delta_G_list[index] = delta_G

        Ni_H2O = volume_fraction_gas/(1-volume_fraction_gas) * (acell_LSCF_slab/2)**3/(kb*T) * p_H2O * 1e5

        a = -1-np.exp(delta_G/(R*T))
        b = Ni_H2O + x + 3
        c = -(x*Ni_H2O + 3*Ni_H2O + 3*x)
        d = 3*x*Ni_H2O
        solutions = cubic_model(a,b,c,d)
        #equation_func = lambda x: equation(x, T=T, p_H2O=  p_H2O, delta_G=delta_G)
        #solution = bisection(0,0.4)

        V_Sr[index] = solutions[0]
    return (V_Sr, delta_G_list)



if __name__ == "__main__":
    data = np.genfromtxt("./lib/sroh2_factsage.csv", delimiter=";")
    T_range = data[100::100,0]
    T = T_range[99]
    mu_H2O = cp_H2O(T, E_DFT_H2O=E_DFT_H2O)
    mu_SrOH2 = cp_SrOH2(T)
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

    fig2.savefig("Ni_H2O_sensitivity.png", dpi=300, transparent=True, format="png")
    plt.show()
