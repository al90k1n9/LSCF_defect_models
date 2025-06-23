"""Functions that are called ubiquitously for many conditions, in all conditions.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy.optimize import curve_fit as cf
import os

local_path = os.path.dirname(os.path.abspath(__file__))
local_path += "/"

default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
N_avagadro = 6.0223*10**23
R = 8.31446261815 #J.K.mol-1

def linear_function(x, a, b):
    return a*x+b

def log_convert(x):
    return np.log(x)

def log_invert(x):
    return np.exp(x)

def delta_oxygen_interpolater(plot=1, pO2 = 0.21):
    """Function to interpolate oxygen understoichometry for a given temperature. The interpolation is valid only within temerature window set by the experimental data. In this case it is the data from Bouwmeester.

    Parameters
    ----------
    plot : bool, optional
        Activate to show plots
    pO2 : float, optional
        Oxygen partial pressure. Defaults to ambient gas composition of 0.21 atm.

    Returns
    -------
    numpy array of size 1x2
        Returns numpy array containing linear interpolation parameters: slope and offset.

    Warnings
    --------
    The function here returns interpolation parameters for the oxygen understoichiometry data used by Federio Monaco for his theis work. The temperature dependance is taken from Bouwmeester. Using this function give you  the parmaeters for a linear function for oxygen understoichometry for temperature dependence that goes to zero at around 930K. This is unphysical. Needs to be changed. Include an assertion error somwhere to avoid extrapolating temperature dependance outside the window defined by Bouwmeester data.

    """
    data873 = np.genfromtxt(local_path + "bouwmeester_873.csv", delimiter=",", skip_header=1)
    data973 = np.genfromtxt(local_path + "bouwmeester_973.csv", delimiter=",", skip_header=1)
    data1073 = np.genfromtxt(local_path + "bouwmeester_1073.csv", delimiter=",", skip_header=1)
    data1173 = np.genfromtxt(local_path + "bouwmeester_1173.csv", delimiter=",", skip_header=1)

    
    popt, pcov = cf(linear_function, data873[:,0], data873[:,1], p0=[0,0])
    popt2, pcov2 = cf(linear_function, data973[:,0], data973[:,1], p0=[0,0])
    popt3, pcov3 = cf(linear_function, data1073[:,0], data1073[:,1], p0=[0,0])
    popt4, pcov4 = cf(linear_function, data1173[:,0], data1173[:,1], p0=[0,0])

    log_pO2 = np.log(pO2)
    T_range = np.asarray([873,973,1073,1173])
    delta_list = np.asarray([linear_function(log_pO2, popt[0], popt[1]), linear_function(log_pO2, popt2[0], popt2[1]), linear_function(log_pO2, popt3[0], popt3[1]), linear_function(log_pO2, popt4[0], popt4[1])])

    interpolated_param, cov = cf(linear_function, T_range, delta_list, p0=[0,0])
    
    max_concentration = 83108 #mol/m**3
    federico_eq_concentration = 82887 #mol/m**3 federico monaco thesis, p.111 tab. 3.7
    federico_delta_eq_concentration = 3*(1-federico_eq_concentration/max_concentration)
    federico_delta_param = [interpolated_param[0], federico_delta_eq_concentration - interpolated_param[0] * 973]
    print("federico monaco delta_oxygen inversion temperature T [K] ", -federico_delta_param[1]/federico_delta_param[0])
    print("federico delta", federico_delta_eq_concentration)

    if plot:
        fig, ax =plt.subplots(layout="constrained")
        ax.plot(data873[:,0], data873[:,1], label="T = 873K", ls="", marker = "s", markerfacecolor="none")
        ax.plot(data973[:,0], data973[:,1], label="T = 973K", ls="", marker = "s", markerfacecolor="none")
        ax.plot(data1073[:,0], data1073[:,1], label="T = 1073K", ls="", marker = "s", markerfacecolor="none")
        ax.plot(data1173[:,0], data1173[:,1], label="T = 1173K", ls="", marker = "s", markerfacecolor="none")
        

        
        artxlist_min = min(data873[0,0], data973[0,0], data1073[0,0], data1173[0,0]) 
        artxlist_max = max(data873[-1,0], data973[-1,0], data1073[-1,0], data1173[-1,0]) 
        artxlist = np.arange(artxlist_min, artxlist_max, (artxlist_max-artxlist_min)/1000)

        secxax = ax.secondary_xaxis("top", functions=(log_invert, log_convert))
        
        ax.plot(artxlist, linear_function(artxlist, popt[0], popt[1]), color=default_colors[0], ls="dashed")
        ax.plot(artxlist, linear_function(artxlist, popt2[0], popt2[1]), color=default_colors[1], ls="dashed")
        ax.plot(artxlist, linear_function(artxlist, popt3[0], popt3[1]), color=default_colors[2], ls="dashed")
        ax.plot(artxlist, linear_function(artxlist, popt4[0], popt4[1]), color=default_colors[3], ls="dashed")

        ax.axvline(x=np.log(0.21), color="black", ls="dotted")
        
        ax.set_xlabel("log pO$_2$")
        ax.set_ylabel("$\\delta$")
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())

        ax.legend(facecolor="none")

        fig2, ax2 = plt.subplots(layout="constrained")
        ax2.plot(T_range, delta_list,  marker="s", ls="", markerfacecolor="none", color="black")
        #artTlist = np.arange(min(T_range), max(T_range), (max(T_range)-min(T_range))/1000)
        artTlist = np.arange(600, 1200)
        ax2.plot(artTlist, linear_function(artTlist, interpolated_param[0], interpolated_param[1]), label="pO2 = " + str(pO2) + " Bouwmeester")
        ax2.plot(artTlist, linear_function(artTlist, federico_delta_param[0], federico_delta_param[1]),  label="pO2 = " + str(pO2) + " Monaco electrode")
        ax2.axvline(x=973, ls="dotted", color="black")
        ax2.axhline(y=0.0227, ls="dotted", color = "black")
        
        ax2.set_xlim(artTlist[0], artTlist[-1])
        ax2.legend()
        ax2.set_xlabel("T [K]")
        ax2.set_ylabel("$\\delta$")
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.yaxis.set_minor_locator(AutoMinorLocator())


        plt.show()
    return federico_delta_param

def pH2_giver(T,p_H2O, p_O2):
    """To estimate hydrogen gas partial pressure assuming equilibrium between water vapour, oxygen gas and hydrogen gas.

    Parameters
    ----------
    T : float
        Temperature in Kelvin
    p_H2O : float
        Water vapour partial pressure
    p_O2 : float
        Oxygen partial presssure.

    Returns
    -------
    float
        Returns hydrogen gas partial pressure

    Raises
    ------
    Uses Shomate equations and therefore Shomate coefficients and the coefficients are given only for range of temperature. If temperature out of this range, raises error!

    """
    #delta mu calculation
    T_0 = 298 #K
    p_0 = 1 #bar


    cp_O2 = 29.4 #K J/mol
    delta_hf298_O2 = 0 #kJ.mol-1
    s0_O2 = 205.2 #J.mol-1.K


    cp_H2 = 28.8 #K J/mol
    delta_hf298_H2 = 0 #kJ.mol-1
    s0_H2O = 130.680 #J.mol-1.K

    cp_H2O = 33.6 #K J/mol
    delta_hf298_H2O = -241.826 #kJ.mol-1
    s0_H2O = 188.835 #J.mol-1.K-1

    t=T/1000
    h_mat = np.array([t,  t**2/2, t**3/3, t**4/4, -1/t, 1, 0, -1])
    s_mat = np.array([np.log(t), t,  t**2/2.0,  t**3/3, -1./(2 * t**2), 0, 1, 0])
    assert ((T<1000 and T>=700) or (T>=1000 and T<1700)), "Temperature out of range in shomate equations"
    if T<1000 and T>=700:
        shomate_coeffs = [[30.09200, 6.832514, 6.793435, -2.534480,0.082139, -250.8810, 223.3967, -241.8264],
                            [33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797, 172.707974, 0.0],
                            [30.03235, 8.772972, -3.988133, 0.788313, -0.741599, -11.32468, 236.1663, 0.0]]
        shomate_coeffs = np.asarray(shomate_coeffs)

    elif T>=1000 and T<1700:
        shomate_coeffs = [[30.09200, 6.832514, 6.793435, -2.534480,0.082139, -250.8810, 223.3967, -241.8264],
                            [18.563083, 12.257357, -2.859786, 0.268238, 1.977990, -1.147438, 156.288133, 0.0],
                            [30.03235, 8.772972, -3.988133, 0.788313, -0.741599, -11.32468, 236.1663, 0.0]]
        shomate_coeffs = np.asarray(shomate_coeffs)
    else:
        return "temperature out of range"
    h = np.dot(shomate_coeffs, h_mat * 1000) #J/mol
    s = np.dot(shomate_coeffs, s_mat) #J K /mol

    stoichiometry = np.asarray([1,-1,-1/2])
    h = h * stoichiometry
    s = s * stoichiometry
    

    hf = np.sum(h) + delta_hf298_H2O * 1000
    gf = hf - T * np.sum(s)

    ph2=p_H2O/np.sqrt(p_O2)*np.exp(gf/(R*T))
    return ph2

def cube_root_complex(z):
    #python has problems taking powers of complex functions. so the necessary function, here
    norm = (z.real**2 + z.imag**2)**(1/6)
    phase = np.arctan(z.imag/z.real)/3
    #print(phase, "phase")
    return complex(norm*np.cos(phase), norm*np.sin(phase))

def cubic_model(a:float,b:float,c:float,d:float, zero_threshold=0):
    """Analytical solution to cubic functions. See wikipedia.

    Parameters
    ----------
    zero_threshold : float, optional
        To accept some leeway in equations with a threshold. Defaults to 0 to perfect analytical solution.
    a : float
        Polynomial coefficient $ax^3 + bx^2 + cx + d = 0$
    b : float
        Polynomial coefficient $ax^3 + bx^2 + cx + d = 0$
    c : float
        Polynomial coefficient $ax^3 + bx^2 + cx + d = 0$
    d : float
        Polynomial coefficient $ax^3 + bx^2 + cx + d = 0$

    Returns
    -------
    numpy array of length 1x3
        Returns an array containing all three solutions


    Warnings
    --------
    Solutions are complex in general. There is no way to order them. You might end up extracting unphysical solution (complex, real outside physical meaning). When strange things are observed, think of this!

    """
    #FUNCTION TO FIND CUBIC ROOTS
    determinant = 18 *a*b*c*d - 4*b**3*d + b**2 * c**2 - 4*a*c**3 - 27*a**2 * d**2
    #print(N, "N")
    if abs(determinant) > zero_threshold: #test for determinant being zero. numerical errors can lead to infinitesimally small values
        delta_zero = b**2 - 3*a*c
        delta_one = 2*b**3 - 9*a*b*c + 27*a**2*d 
        inner_sqrt = delta_one**2 - 4 * delta_zero**3
        re_inner_sqrt = inner_sqrt
        if inner_sqrt >= 0:
            c_minus = np.cbrt((delta_one - np.sqrt(delta_one**2 - 4 * delta_zero**3))/2)
            c_plus = np.cbrt((delta_one + np.sqrt(delta_one**2 - 4 * delta_zero**3))/2)
        else:
            inner_sqrt = np.abs(inner_sqrt)
            c_minus = cube_root_complex(delta_one/2 - complex(0,np.sqrt(inner_sqrt)/2))
            c_plus = cube_root_complex(delta_one/2 - complex(0,np.sqrt(inner_sqrt)/2))


        xi = complex(-1/2, np.sqrt(3)/2)
        xi_coeffs = np.asarray([1,xi, 1/xi])
        xi_powers = np.asarray([xi_coeffs**0, xi_coeffs, xi_coeffs**2])
        terms_minus = -1/(3*a)*np.asarray([b, c_minus, delta_zero/c_minus])
        terms_plus = -1/(3*a)*np.asarray([b, c_plus, delta_zero/c_plus])
        solutions_plus =np.dot( xi_powers, terms_plus)
        solutions_minus = np.dot(xi_powers, terms_minus)
        #print(solutions_minus[0] -solutions_plus[0], "difference in solutions ", solutions_minus[0])
        return solutions_plus
    else:
        D0 = c**2 -3*b*d
        D1 = 2*c**3 - 9 *b*c*d + 27*a*d**2
        if abs(D0)<zero_threshold and abs(D1)<zero_threshold:
            return np.asarray([-b/(3*a),-b/(3*a),-b/(3*a)])
        else:
            # Double root case
            double_root = -b/(3*a) + D1/D0
            # Find the single root using the fact that sum of roots = -b/a
            single_root = -b/a - 2*double_root
            return np.array([double_root, double_root, single_root])


def quadratic_model(a:float,b:float,c:float,x0=0.4):
    """Analytical solution of quadratic polynomial

    Parameters
    ----------
    x0 : float, optional
        Initial amounts of Sr in LSCF in mol per unit LSCF volume. Defaults to 0.4. Used to select the physical solution, postive and below max amount of Sr.
    a : float
        Polynomial coefficient $ax^2 + bx + c = 0$
    b : float
        Polynomial coefficient $ax^2 + bx + c = 0$
    c : float
        Polynomial coefficient $ax^2 + bx + c = 0$

    Returns
    -------
    float
        Returns the one physical solutions.

    Warnings
    --------
    When both solutions are unphysical, it returns 2. This should be visible in graphs and return here in that case. Create an error otherwise maybe.

    """
    delta = b**2 - 4*c*a
    solution1 = (-b + np.sqrt(delta))/(2*a)
    solution2 = (-b - np.sqrt(delta))/(2*a)
    if solution1<=x0 and solution1>=0:
        return solution1
    elif solution2 <=x0 and solution2>=0:
        return solution2
    else:
        print(solution1, solution2)
        return 2

def bisection(xmin:float, xmax:float, func, maxiter=10000, tol=1e-20):
   """Bisection method for numerical resolution of an equation

   Parameters
   ----------
   xmin : float
       Lower bound of the window in which the solution is present.
   xmax : float
       Upper bouund of the window in which the solution is present.
   func : function
       Function for which the solution needs to be found.
   maxiter : int
       Maximum iteration to not exceed, to avoid infinite loops.
   tol : float
       Precision of the solution.

   Returns
   -------
   float
       Returns the solution.J

   """
   fmin = func(xmin)
   fmax = func(xmax)
   if  fmin * fmax > 0:
       raise ValueError("no solution in the given interval")
   for iteration in range(maxiter):
       c= (xmin+xmax)/2
       if abs(func(c))<tol or (xmax-xmin)/2<tol:
           return round(c, 7)
       
       if func(xmin)*func(c)<0:
           xmax = c
       else:
           xmin = c
   raise ValueError("maxiter reached")

from lib.dft_energies_0K import E_DFT_H2O, ev2J_p_mol
from lib.chemical_potentials import chem_pot_H2O   
kb = 1.380649 * 10**(-23) #J/K
m_H2O = 18.01528 / (N_avagadro*1000) #in kg
hbar = 1.054571817*10**(-34) #reduced planck's constant in J.s


def surface_coverage_H2O(T, x_H2O, E_ads, P, chem_pot = 0):
    """Estimate surface coverage of adsborbed water molecules.

    Parameters
    ----------
    T : float, numpy array
        Temperature or temperature window in Kelvins
    x_H2O : float
        Molar fraction of water vapour.
    E_ads : float
        Adsorption energy water molecule onto LSCF surface in J/mol.
    P : float
        Total pressure of the system.
    chem_pot : bool
        Boolean to determine which water vapour chemical potential to use. Defaults to using experimental water vapour chemical potential. When passed 1, chemical potential with only translational degrees of freedom is considered.

    Returns
    -------
    tuple of length 2
        Returns a tuple containing two elements: surface coverage in the same type as passed T, water vapour chemical potential in J/mol used for the calculation.

    Warnings
    --------
    If using experimental chemical potential, instead of a model chemical potential that considers only the translational degrees of freedom of an ideal gas, make sure to include configurational entropy of water vapour in the Gibbs free energy term, by passing an argument to the pressure in the chemical potential term, instead of keeping it at the default value of 1 atm (reference state).

    """
    #chem_pot 0 for experimental chemical potential
    #chem_pot 1 for translational chemical potential
    if chem_pot == 0: 
        chemical_potential = chem_pot_H2O(T, E_DFT_H2O=0, P=x_H2O*P) #J/mol
    else:
        chemical_potential = N_avagadro * kb*T* np.log(x_H2O*P*1e5/(kb*T) * (2*np.pi*hbar**2/(m_H2O*kb*T))**(3/2)) #J/mol
    #BE CAREFUL YOU NEED TO GET THE CHEMICAL POTENTIAL CORRECTION; IT'S IMPORTANT TO PUT THE DFT ENERGY PARAMETER TO ZERO
    exp_term = np.exp(-(-E_ads+chemical_potential)/(R*T))
    theta = 1/(1+exp_term)
    return (theta, chemical_potential)
