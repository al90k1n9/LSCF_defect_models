"""SrO formation models in humid air conditions, using hydrogen gas.
All of the cases presented in this file correspond to Sr coming from the LSCF surface, unless mentioned otherwise.
See lib.SrO_hydroxylated_models for cases used adsorbed water molecules.
"""

import sys
sys.path.append("H:\Documents\SrO_defect_model")
from lib.chemical_potentials import *
from lib.dft_energies_0K import * #importing all the values of dft energies
from lib.auxilliary_functions import *
from H_vibration import *

def case1(T_range, x0=0.4, x_O2 = 0.21, x_H2O = 0.08, P=1, sensitivity_shift = 0):
    """Case where the remanining hydrogen atoms are recombined to form hydrogen gas.

    Parameters
    ----------
    T_range : numpy array
        Temperature window in Kelvins
    x0 : float, optional
        Initial amounts of Sr in moles per unit LSCF volume. Defaults to standard LSCF composition 0.4.
    x_O2 : float, optional
        Molar fraction of oxygen gas. Defaults to ambient gas composition of 0.21.
    x_H2O : float, optional
        Water vapour molar fraction. Defaults to 0.08, the experimental reference used by Sasson et al.
    P : float, optional
        Total pressure of the system. Defaults to ambient gas total pressure of 1 atm/Bar.

    Returns
    -------
    typle of length 3
        Tuple of three elements: all numpy arrays : Sr vancancies in moles per unit LSCF volume, delta_G in J/mol, hydrogen gas partial pressure in atm as a function of temperature.

    """
    p_O2 = x_O2 * P
    p_H2O = x_H2O * P
    V_Sr= []
    delta_E = (E_LSCF_slab_Sr_vac_surf + 2*E_SrO_epitax + 2*E_DFT_H2- (E_LSCF_slab + 2*E_DFT_H2O))/2 + E_int + sensitivity_shift
    #print(delta_E/ev2J_p_mol, " delta_E of case 1")
    delta_G_range = []
    p_H2_list = []
    for T in T_range:
        p_H2 = pH2_giver(T, p_H2O, p_O2)
        p_H2_list.append(p_H2)

        delta_G = (E_LSCF_slab_Sr_vac_surf + 2*chem_pot_SrO(T) + 2*chem_pot_H2(T, E_DFT_H2, P=P)- (E_LSCF_slab + 2*chem_pot_H2O(T, E_DFT_H2O, P=P)))/2 + E_int + sensitivity_shift
        delta_G_range.append(delta_G)
        #if T == 1000: print("delta G at T = 1000", delta_G/ev2J_p_mol)
        K = np.exp(-delta_G/(R*T))

        #DEFINITION OF THIRD DEGREE POLYNOMIAL THAT IS REACTION MECHANISM SPECIFIC
        #SHOULD BE VERIFIED EVERY TIME THE REACTION MECHANISM IS CHANGED
        N = K * p_H2O/p_H2 #notice that the total pressure cancels out in the case
        a = 4+4*N
        b= 4*(x0+1*N)
        c= x0**2 + (1-x0)*N * (1+3*x0)
        d = -N * x0 * (1-x0)**22 
        solution= cubic_model(a,b,c,d)
        #print(x0_minus- x0_plus)
        #print(equation(x0_minus), equation(x0_plus))
        V_Sr.append(solution[0])
        #print(solution[0])
    return (np.asarray(V_Sr), np.asarray(delta_G_range), np.asarray(p_H2_list))

def case2(T_range, x0=0.4, x_H2O = 0.08, P=1,sensitivity_shift = 0):
    """Case where the remaining two hydrogen are stablised in the Sr vacancy created by SrO formation.

    Parameters
    ----------
    T_range : numpy array
        Temperature window in Kelvins.
    x0 : float, optional
        Initial amounts of Sr in mol per unit LSCF volume. Defaults to standard LSCF composition of 0.4.
    x_H2O : float, optional
        Water vapour molar fraction. Defaults to 0.08, value used experimentally by Sassone et al.
    P : float, optional
        Total pressure in the gas. Defaults to ambient total pressure of 1 atm/Bar.

    Returns
    -------
    tuple of length 2
        Returns a tuple containing 2 numpy array: Sr vacancies in mol per unit LSCF volume, delta_G n J/mol as a function of temperature given in T_range.

    """
    p_H2O = x_H2O * P
    V_Sr= []
    delta_E = (2*E_SrO_epitax + E_LSCF_double_hydrogenated - (E_LSCF_slab+ 2 * E_DFT_H2O)) / 2 + E_int + sensitivity_shift
    #print(delta_E/ev2J_p_mol, " delta_E of case 2")
    delta_G_range = []
    for T in T_range:
        delta_G = (2*chem_pot_SrO(T) + E_LSCF_double_hydrogenated - (E_LSCF_slab+ 2 *chem_pot_H2O(T, E_DFT_H2O, P=P))) / 2  + E_int + 2*OH_bond_vibration + sensitivity_shift
        delta_G_range.append(delta_G)
        #if T == 1000: print("case 2: delta G at T = 1000", delta_G/ev2J_p_mol)
        K = np.exp(-delta_G/(R*T)) 
        N = K * p_H2O/P 
        V_Sr.append(N/(1+N)*x0)

    return (np.asarray(V_Sr), np.asarray(delta_G_range))

def case3(T_range, x0=0.4, x_O2 = 0.21, x_H2O = 0.08, P=1):
    """Case where one of the hydrogen in stabilised in the Sr vacancy left by Sr formation and other hydrogen atom forms hydrogen gas.

    Parameters
    ----------
    T_range : numpy array
        Temperature window in Kelvins.
    x0 : float, optional
        Initial amounts of Sr in moles per unit LSCF volume. Defaults to standard LSCF composition of 0.4.
    x_O2 : float, optional
        Oxygen  gas molar fraction. Defaults to ambient gas composition of 0.21.
    x_H2O : float, optional
        Water vapour partial pressure. Defaults to 0.08, one used by Sassone et al.
    P : float, optional
        Total pressure. Defaults to 1 atm/Bar, ambient gas total pressure.

    Returns
    -------
    Tuple of length 2
        Returns a tuple containing 2 numpy arrays: Sr vacancies in moles per unit LSCF volume, delta_G in J/mol as a function of temperature.

    """
    p_O2 = x_O2 * P
    p_H2O = x_H2O * P
    V_Sr= []
    delta_E = (E_LSCF_single_hydrogenated + 2*E_SrO_epitax + E_DFT_H2- (E_LSCF_slab + 2*E_DFT_H2O))/2 + E_int
    print(delta_E/ev2J_p_mol, " delta_E of case 3")
    delta_G_range = []
    for T in T_range:
        p_H2 = pH2_giver(T, p_H2O, p_O2)

        delta_G = (E_LSCF_single_hydrogenated + 2*chem_pot_SrO(T) + chem_pot_H2(T, E_DFT_H2, P=P) - (E_LSCF_slab + 2*chem_pot_H2O(T, E_DFT_H2O, P=P)))/2 + E_int
        delta_G_range.append(delta_G)
        if T == 1000: print("case 3: delta G at T = 1000", delta_G/ev2J_p_mol)
        K = np.exp(-delta_G/(R*T))

        #DEFINITION OF THIRD DEGREE POLYNOMIAL THAT IS REACTION MECHANISM SPECIFIC
        #SHOULD BE VERIFIED EVERY TIME THE REACTION MECHANISM IS CHANGED
        N = K * p_H2O/P * np.sqrt(P/p_H2)
        a = 1-N
        b = x0+N
        c = -N*(x0-x0**2)
        V_Sr.append(quadratic_model(a,b,c))
    return (np.array(V_Sr), np.asarray(delta_G_range))

def case4(T_range, x0=0.4, x_O2 = 0.21, x_H2O = 0.08, P=1):
    """Same as case 1 but with Sr coming from the bulk of LSCF.

    Parameters
    ----------
    T_range : numpy array
        Temperature window in Kelvins
    x0 : float, optional
        Initial amounts of Sr in moles per unit LSCF volume. Defaults to standard LSCF composition 0.4.
    x_O2 : float, optional
        Molar fraction of oxygen gas. Defaults to ambient gas composition of 0.21.
    x_H2O : float, optional
        Water vapour molar fraction. Defaults to 0.08, the experimental reference used by Sasson et al.
    P : float, optional
        Total pressure of the system. Defaults to ambient gas total pressure of 1 atm/Bar.

    Returns
    -------
    typle of length 3
        Tuple of three elements: all numpy arrays : Sr vancancies in moles per unit LSCF volume, delta_G in J/mol, hydrogen gas partial pressure in atm as a function of temperature.

    """
    p_O2 = x_O2 * P
    p_H2O = x_H2O * P
    V_Sr= []
    delta_E = (E_LSCF_slab_Sr_vac_bulk + E_DFT_H2 + E_SrO_epitax) - (E_LSCF_slab + E_DFT_H2O) + E_int
    print(delta_E/ev2J_p_mol, " delta_E of case 4")
    delta_G_range = []
    for T in T_range:
        p_H2 = pH2_giver(T, p_H2O, p_O2)

        delta_G = (E_LSCF_slab_Sr_vac_bulk + chem_pot_H2(T, E_DFT_H2, P=P) + chem_pot_SrO(T)) - (E_LSCF_slab + chem_pot_H2O(T, E_DFT_H2O, P=P)) + E_int
        delta_G_range.append(delta_G)
        if T == 1000: print("delta G at T = 1000", delta_G/ev2J_p_mol)
        K = np.exp(-delta_G/(R*T))

        #DEFINITION OF THIRD DEGREE POLYNOMIAL THAT IS REACTION MECHANISM SPECIFIC
        #SHOULD BE VERIFIED EVERY TIME THE REACTION MECHANISM IS CHANGED
        N = K * p_H2O/p_H2 #notice that the total pressure cancels out in the case
        a = 4+4*N
        b= 4*(x0-1*N)
        c= x0**2 + (1-x0)*N * (1+3*x0)
        d = -N * x0 * (1-x0)**2
        solution= cubic_model(a,b,c,d)
        #print(x0_minus- x0_plus)
        #print(equation(x0_minus), equation(x0_plus))
        V_Sr.append(solution[0])
        #print(solution[0])
    return (V_Sr, delta_G_range)


def case5(T_range, x0 = 0.4, x_H2O = 0.08, P = 1):
    """Same as case 1 but here the hydrogen partial pressure is deduced by assuming a fixed volume of gas, that is the pores of LSCF.

    Parameters
    ----------
    T_range : numpy array
        Temperature window in Kelvins
    x0 : float, optional
        Initial amounts of Sr in moles per unit LSCF volume. Defaults to 0.4, the standard LSCF composition
    x_H2O : float, optional
        Water vapour molar fraction. Defaults to 0.08, value used experimentally in Sassone et al.
    P : float, optional
        Total pressure. Defaults to ambient gas total pressure of 1 atm/bar.

    Returns
    -------
    tuple of length 2
        Tuple containing 2 numpy arrays: Sr vacancies in moles per unit LSCF volume, delta_G in J/mol as a function of temperature.

    """
    #case similar to 1 in KV notation. but the activity is expressed differently
    #experimentally equivalent to keep a lid on the LSCF surface and using only the porosity of the LSCF as gas volume
    V_Sr= []
    delta_G_list = []
    delta_E = (E_LSCF_slab_Sr_vac_surf + 2*E_SrO_epitax + 2*E_DFT_H2- (E_LSCF_slab + 2*E_DFT_H2O))/2 + E_int
    print("delta E ", delta_E/ev2J_p_mol)
    for T in T_range:
        delta_G = (E_LSCF_slab_Sr_vac_surf + 2*chem_pot_SrO(T) + 2*chem_pot_H2(T, E_DFT_H2, P=P)- (E_LSCF_slab + 2*chem_pot_H2O(T, E_DFT_H2O, P=P)))/2 + E_int
        delta_G_list.append(delta_G)
    return (np.asarray(V_Sr), np.asarray(delta_G_list))


def ph2_sensitivity_case1(x_H2_range, T=1000, x0 = 0.4, x_H2O=0.08, P=1):
    """Sensitivity analysis to see the effect of hydrogen partial pressure in case 1.

    Parameters
    ----------
    x_H2_range : numpy array
        Hydrogen gas molar fraction array.
    T : float, optional
        Temperature in Kelvin. Defaults to 1000K.
    x0 : float, optional
        Initial amounts of Sr in mol per unit LSCF volume. Defaults to standard LSCF composition of 0.4.
    x_H2O : float, optional
        Water vapour molar fraction. Defaults to 0.08, value used experimentally by Sassone et al.
    P : float, optional
        Total pressures. Defaults to 1 atm/Bar as in ambient gas total pressure.

    Returns
    -------
    numpy array
        Returns Sr vacancies in moles per unit LSCF volume

    """
    p_H2O = x_H2O * P
    V_Sr= []
    delta_G = (E_LSCF_slab_Sr_vac_surf + 2*chem_pot_SrO(T) + 2*chem_pot_H2(T, E_DFT_H2, P=P) - (E_LSCF_slab + 2*chem_pot_H2O(T, E_DFT_H2O, P=P)))/2 + E_int
    K = np.exp(-delta_G/(R*T))
    for x_H2 in x_H2_range:
        p_H2 = x_H2 * P
        N = K * p_H2O/p_H2 #notice that the total pressure cancels out in the case
        a = 4+4*N
        b= 4*(x0-1*N)
        c= x0**2 + (1-x0)*N * (1+3*x0)
        d = -N * x0 * (1-x0)**2
        solution= cubic_model(a,b,c,d)
        V_Sr.append(solution[0])
    return np.asarray(V_Sr)
