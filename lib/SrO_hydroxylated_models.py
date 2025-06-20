"""SrO formation in humid conditions, now with starting hydroxylated surfaces i.e. water adsorbed onto LSCF surfaces.
"""

from lib.chemical_potentials import *
from lib.dft_energies_0K import * #importing all the values of dft energies
from lib.auxilliary_functions import *



def case1(T_range, x=0.4, x_O2 = 0.21, x_H2O = 0.08, P=1):
    """Case where the remaining two hydrogen atoms are recombined to form hydrogen gas.

    Parameters
    ----------
    T_range : numpy array
        Temperature window in Kelvins
    x : float, optional
        Initial amounts of Sr in mols per unit LSCF volume. Defaults to 0.4.
    x_O2 : float, optional
        Oxygen gas molar fraction. Defaults to ambient gas composition of 0.21.
    x_H2O : float, optional
        Water vapour molar fraction. Defaults to 0.08 that is used experimentally by Sassone et al.
    P : float, optional
        Total pressure. Defaults to 1 atm/Bar that is ambient gas total pressure.

    Returns
    -------
    tuple of length 3
        Returns a tuple containing 3 numpy arrays: Sr vacancies in mols per unit LSCF volume, delta_G in J/mol, surface coverage as a function of temperature given in T_range.

    """
    p_H2O = x_H2O * P
    p_O2 = x_O2 * P
    V_Sr= []
    delta_E = (E_LSCF_slab_Sr_vac_surf + 2*E_SrO_epitax + 2*E_DFT_H2- E_LSCF_hydroxilated )/2 + E_int
    #print(delta_E/ev2J_p_mol, " of case 1")
    delta_G_list = []
    theta_list = []
    for T in T_range:
        p_H2 = pH2_giver(T, p_H2O, p_O2)

        delta_G = (E_LSCF_slab_Sr_vac_surf + 2*chem_pot_SrO(T) + 2*chem_pot_H2(T, E_DFT_H2, P=P)- E_LSCF_hydroxilated )/2 + E_int
        delta_G_list.append(delta_G)
        theta = surface_coverage_H2O(T,p_H2O/P, E_ads, P)[0]
        N = theta/(1-theta) * np.exp(-delta_G/(R*T)) * P/p_H2

        a = 4+4*N
        b= 4*(x-1*N)
        c= x**2 + (1-x)*N * (1+3*x)
        d = -N * x * (1-x)**2
        solution= cubic_model(a,b,c,d)
        V_Sr.append(solution[0])
        theta_list.append(theta)
    return (np.asarray(V_Sr), np.asarray(delta_G_list), np.asarray(theta_list))

def case2(T_range, x=0.4, x_H2O = 0.08, P=1, sensitivity_shift = 0):
    """Case where the two remaining hydrogen atoms are stabilised in Sr vacancy left by SrO formation.

    Parameters
    ----------
    T_range : numpy array
        Temperature window in Kelvins
    x : float, optional
        Initial amounts of Sr in mols per unit LSCF volume. Defaults to 0.4, the standard composition of LSCF.
    x_H2O : float, optional
        Water vapour molar fraction. Defaults to 0.08, value experimentally used by Sassone et al
    P : float, optional
        Total pressure. Defaults to ambient gas total pressure 1 atm/Bar.
    sensitivity_shift : float
        Arbitrary shift in J/mol introduced in delta_G to do sensitivity analysis.

    Returns
    -------
    tuple of length 3
        Returns a tuple containing 3 numpy arrays: Sr vacancies in mols per unit LSCF volume, delta_G in J/mol, surface coverage of water molecules as a function of temperature.

    Warnings
    --------
    Keep the sensitivity_shift value to its default value of 0 for physically meaningful results.

    """
    p_H2O = x_H2O * P
    V_Sr= []
    delta_E = (E_LSCF_double_hydrogenated + 2*E_SrO_epitax - E_LSCF_hydroxilated )/2 + E_int + sensitivity_shift
    #print(delta_E/ev2J_p_mol, " of case 2")
    
    delta_G_list = []
    theta_list = []
    for T in T_range:
        delta_G = (E_LSCF_double_hydrogenated + 2*chem_pot_SrO(T) - E_LSCF_hydroxilated )/2 + E_int + sensitivity_shift
        delta_G_list.append(delta_G)
        theta, chemical_potential_used = surface_coverage_H2O(T,p_H2O/P, E_ads, P)
        N = theta/(1-theta) * np.exp(-delta_G/(R*T))
        V_Sr.append(N/(1+N)*x)
        theta_list.append(theta)
    return (np.asarray(V_Sr), np.asarray(delta_G_list), np.asarray(theta_list))

def case3(T_range, x=0.4, x_O2 = 0.21, x_H2O = 0.08, P=1):
    """One of the reamining hydrogen atoms is stabilised in the Sr vacancy and the other one is used to produce hydrogen gas.

    Parameters
    ----------
    T_range : numpy array
        Temperature window in Kelvins
    x : float, optional
        Initial amounts of Sr in mols per unit LSCF volume. Defaults to standard LSCF composition of 0.4.
    x_O2 : float, optional
        Oxygen gas molar fraction. Defaults to ambient gas composition of 0.21.
    x_H2O : float, optional
        Water vapour molar fracation. Defaults to 0.08, experimentally used value by Sassone et al.
    P : float, optional
        Total pressure. Defaults to 1, ambient gas total pressure.

    Returns
    -------
    tuple of length 3
        Returns a tuple containing 3 numpy arrays: Sr vacacnies in mols per unit LSCF volume, delta_G in J/mol, surface coverage as a function of temperature given in T_range.

    """
    p_O2 = x_O2 * P
    p_H2O = x_H2O * P
    V_Sr= []
    delta_E = (E_LSCF_single_hydrogenated + 2*E_SrO_epitax + E_DFT_H2 - E_LSCF_hydroxilated )/2 + E_int
    #print(delta_E/ev2J_p_mol, " of case 3")
    delta_G_list = []
    theta_list= []
    for T in T_range:
        p_H2 = pH2_giver(T, p_H2O, p_O2)

        delta_G = (E_LSCF_single_hydrogenated + 2*chem_pot_SrO(T) + chem_pot_H2(T, E_DFT_H2, P=P) - E_LSCF_hydroxilated )/2 + E_int
        delta_G_list.append(delta_G)

        theta = surface_coverage_H2O(T,p_H2O/P, E_ads, P)
        N = theta/(1-theta) * np.exp(-delta_G/(R*T)) * np.sqrt(P/p_H2)

        a = 1-N
        b = x+N
        c = -N*(x-x**2)
        V_Sr.append(quadratic_model(a,b,c))
        theta_list.append(theta)
    return (V_Sr, np.asarray(delta_G_list), np.asarray(theta_list))



def ph2_sensitivity_case1(x_H2_range, T=1000, x0 = 0.4, x_H2O=0.08, P=1):
    """Enables the sensitivity analysis of hydrgen partial pressure on case1

    Parameters
    ----------
    x_H2_range : numpy array
        Range of hydrogen gas molar fractions.
    T : float, optional
        Temperature in Kelvin. Defaults to 1000K
    x0 : float, optional
        Initial amounts of Sr in mols per unit LSCF volume. Defaults to standard LSCF composition of 0.4.
    x_H2O : float, optional
        Water vapour molar fraction. Defaults to 0.08, experimental value used by Sassone et al.
    P : float, optional
        Total pressures. Defaults to ambient gas total pressure of 1 atm/Bar.

    Returns
    -------
    numpy array
        Sr vacancies in mols per unit LSCF volume as a function of temperature provided in T_range.

    """
    p_H2O = x_H2O * P
    V_Sr= []
    delta_G = (E_LSCF_slab_Sr_vac_surf + 2*chem_pot_SrO(T) + 2*chem_pot_H2(T, E_DFT_H2, P=P)- E_LSCF_hydroxilated )/2 + E_int
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

