"""Different cases of SrO formation in dry air conditions
"""

from lib.chemical_potentials import *
from lib.dft_energies_0K import * #importing all the values of dft energies
from lib.auxilliary_functions import *

def case1(T_range, x=0.4, x_O2 = 0.21, P=1, sensitivity_shift = 0):
    """Case with Sr from surface and oxygen from the atmosphere

    Parameters
    ----------
    T_range : numpy array
        Temperature window in Kelvins.
    x : float, optional
        Initial amount of Sr in mol per unit LSCF volume in the material. Defaults to 0.4 standard LSCF stoichiometry
    x_O2 : float, optional
        Oxygen gas molar fraction. Defaults to ambient gas composition of 0.21
    P : float, optional
        Total pressure. Defaults to ambient gas total pressure of 1 atm/Bar.
    sensitivity_shift : float, optional
        Arbitrary shift in energy to add in delta_G for sensitivity analysis. Defaults to zero for meaninful results.

    Returns
    -------
    tuple of length 2
        Returns two numpy arrays: Sr vacancies in unit LSCF volume, reference Gibbs free energy in J/mol as a function of temperature given in T_range

    Warnings
    --------
    Sensitivity analysis argument should be kept to its default value for physically relvant results.

    """
    p_O2 = P*x_O2
    V_Sr= []
    delta_G_list = []
    delta_E = (E_LSCF_slab_Sr_vac_surf + 2 * E_SrO_epitax - (E_LSCF_slab + E_DFT_O2))/2 + E_int + sensitivity_shift
    for T in T_range:
        delta_G = (E_LSCF_slab_Sr_vac_surf + 2 * (chem_pot_SrO(T)) - (E_LSCF_slab + chem_pot_O2(T, E_DFT_O2, P=P)))/2 + E_int + sensitivity_shift
        #if T==973: print("case 1", delta_G/ev2J_p_mol)
        K = np.exp(-delta_G/(R*T))
        #print(K)

        #DEFINITION OF THIRD DEGREE POLYNOMIAL THAT IS REACTION MECHANISM SPECIFIC
        #SHOULD BE VERIFIED EVERY TIME THE REACTION MECHANISM IS CHANGED
        N = K * np.sqrt(p_O2/P)
        a = 4+4*N
        b= 4*(x-1*N)
        c= x**2 + (1-x)*N * (1+3*x)
        d = -N * x * (1-x)**2
        solution= cubic_model(a,b,c,d)
        #print(x0_minus- x0_plus)
        #print(equation(x0_minus), equation(x0_plus))
        V_Sr.append(solution[0])
        delta_G_list.append(delta_G)
        #print(solution[0])
    return (np.asarray(V_Sr),np.asarray(delta_G_list))

def case2(T_range, x=0.4, x_O2 = 0.21, P=1, sensitivity_shift = 0):
    """Case with Sr coming from LSCF bulk and oxygen from atmosphere.

    Parameters
    ----------
    T_range : numpy array
        Temperature window in Kelvins.
    x : float, optional
        Initial amount of Sr in mol per unit LSCF volume in the material. Defaults to 0.4 standard LSCF stoichiometry
    x_O2 : float, optional
        Oxygen gas molar fraction. Defaults to ambient gas composition of 0.21
    P : float, optional
        Total pressure. Defaults to ambient gas total pressure of 1 atm/Bar.
    sensitivity_shift : float, optional
        Arbitrary shift in energy to add in delta_G for sensitivity analysis. Defaults to zero for meaninful results.

    Returns
    -------
    tuple of length 2
        Returns two numpy arrays: Sr vacancies in unit LSCF volume, reference Gibbs free energy in J/mol as a function of temperature given in T_range

    Warnings
    --------
    Sensitivity analysis argument should be kept to its default value for physically relvant results.

    """
    p_O2 = x_O2*P
    delta_G_list =[]
    V_Sr= []
    for T in T_range:
        delta_G = E_LSCF_slab_Sr_vac_bulk + chem_pot_SrO(T) - (E_LSCF_slab + 0.5*(chem_pot_O2(T, E_DFT_O2, P=P))) + E_int + sensitivity_shift
        if T==973: print("case 2", delta_G/ev2J_p_mol)
        delta_G_list.append(delta_G)
        K = np.exp(-delta_G/(R*T))
        #print(K)

        #DEFINITION OF THIRD DEGREE POLYNOMIAL THAT IS REACTION MECHANISM SPECIFIC
        #SHOULD BE VERIFIED EVERY TIME THE REACTION MECHANISM IS CHANGED
        N = K * np.sqrt(p_O2/P)
        a = 4+4*N
        b= 4*(x-1*N)
        c= x**2 + (1-x)*N * (1+3*x)
        d = -N * x * (1-x)**2
        solution= cubic_model(a,b,c,d)
        V_Sr.append(solution[0])
    return (np.asarray(V_Sr),np.asarray(delta_G_list))

def case3(T_range, x=0.4, delta_oxygen_parameters=np.asarray([0,0]), sensitivity_shift = 0):
    """Case with Sr from surface and oxygen from LSCF sub surface.

    Parameters
    ----------
    T_range : numpy array
        Temperature window in Kelvins
    x : float, optional
        Initial amount of Sr in mol per unit LSCF volume in LSCF. Defaults to 0.4
    delta_oxygen_parameters : numpy array for size 1x2, optional
        Linear interpolation parameters; slope and the intercept, to calculate the oxygen understoichometry as a function of any given temperature, for a given oxygen partial pressure.
    sensitivity_shift : float, optional
        Arbitrary shift in energy to add in delta_G for sensitivity analysis. Defaults to zero for meaninful results.


    Returns
    -------
    tuple of length 3
        Returns 3 numpy arrays: Sr vacancies in mol per unit LSCF volume, delta_G in J/mol, oxygen understoichiometry in mol per unit LSCF volume as a function of temperature given in T_range.

    Warnings
    --------
    Sensitivity analysis argument should be kept to its default value for physically relvant results.

    """
    delta_G_list = []
    delta_oxygen_list = []
    V_Sr=[]
    delta_E = E_LSCF_slab_Sr_surf_O_sub_surf /2 +  E_SrO_epitax - (E_LSCF_slab/2) + E_int + sensitivity_shift
    delta_oxygen_parameter_condition = delta_oxygen_parameters[0] != 0 and delta_oxygen_parameters[1] != 0
    if delta_oxygen_parameter_condition: inversion_temperature = -delta_oxygen_parameters[1]/delta_oxygen_parameters[0]
    for T in T_range:

        delta_G = E_LSCF_slab_Sr_surf_O_sub_surf /2 +  chem_pot_SrO(T) - (E_LSCF_slab/2) + E_int + sensitivity_shift
        if T==973: print("case 3", delta_G/ev2J_p_mol)
        delta_G_list.append(delta_G)
        #delta_oxygen = delta_oxygen_parameters[0] * T + delta_oxygen_parameters[1]
        #if delta_oxygen<0: delta_oxygen = 0 #understoichiometry cannot be negative.
        if delta_oxygen_parameter_condition: delta_oxygen = delta_oxygen_parameters[0] * np.log(1+np.exp(T-inversion_temperature))
        else: delta_oxygen = 0
        delta_oxygen_list.append(delta_oxygen)
        K = np.exp(-delta_G/(R*T))
        a = (1-(1/K))
        b = -(3-delta_oxygen)-x - delta_oxygen/K
        c = (3-delta_oxygen)*x

        V_Sr.append(quadratic_model(a,b,c,x))
    return(np.asarray(V_Sr), np.asarray(delta_G_list), np.asarray(delta_oxygen_list))

def case4(T_range, x=0.4, delta_oxygen_parameters = [0,0], sensitivity_shift = 0):
    """Case with Sr coming from LSCF bulk and oxygen from LSCF sub surface.

    Parameters
    ----------
    T_range : numpy array
        Temperature window in Kelvins
    x : float, optional
        Initial amount of Sr in mol per unit LSCF volume in LSCF. Defaults to 0.4
    delta_oxygen_parameters : numpy array for size 1x2, optional
        Linear interpolation parameters; slope and the intercept, to calculate the oxygen understoichometry as a function of any given temperature, for a given oxygen partial pressure.
    sensitivity_shift : float, optional
        Arbitrary shift in energy to add in delta_G for sensitivity analysis. Defaults to zero for meaninful results.


    Returns
    -------
    tuple of length 3
        Returns 3 numpy arrays: Sr vacancies in mol per unit LSCF volume, delta_G in J/mol, oxygen understoichiometry in mol per unit LSCF volume as a function of temperature given in T_range.

    Warnings
    --------
    Sensitivity analysis argument should be kept to its default value for physically relvant results.

    """
    delta_G_list = []
    V_Sr=[]
    delta_E = E_LSCF_slab_SrO_bulk +  E_SrO_epitax - (E_LSCF_slab) + E_int + sensitivity_shift
    delta_oxygen_list = []
    delta_oxygen_parameter_condition = delta_oxygen_parameters[0] != 0 and delta_oxygen_parameters[1] != 0
    if delta_oxygen_parameter_condition: inversion_temperature = -delta_oxygen_parameters[1]/delta_oxygen_parameters[0]
    for T in T_range:
        delta_G = E_LSCF_slab_SrO_bulk +  chem_pot_SrO(T) - (E_LSCF_slab) + E_int + sensitivity_shift
        #if T==973: print("case 4", delta_G/ev2J_p_mol)
        delta_G_list.append(delta_G)
        #delta_oxygen = delta_oxygen_parameters[0] * T + delta_oxygen_parameters[1]
        #if delta_oxygen<0: delta_oxygen = 0 #understoichiometry cannot be negative.
        if delta_oxygen_parameter_condition: delta_oxygen = delta_oxygen_parameters[0] * np.log(1+np.exp(T-inversion_temperature))
        else: delta_oxygen = 0
        delta_oxygen_list.append(delta_oxygen)
        K = np.exp(-delta_G/(R*T))
        a = (1-(1/K))
        b = -(3-delta_oxygen)-x - delta_oxygen/K
        c = (3-delta_oxygen)*x

        V_Sr.append(quadratic_model(a,b,c,x))
    return(np.asarray(V_Sr), np.asarray(delta_G_list), np.asarray(delta_oxygen_list))

def case5(T_range, x=0.4, x_O2 = 0.21, P=1):
    """Case where LSCF bulk systems are used instead of slabs. Oxygen comes from atmosphere.

    Parameters
    ----------
    T_range : numpy array
        Temperature window in Kelvins
    x : float, optional
        Initial amounts of Sr in mol per unit LSCF volume. Defaults to 0.4
    x_O2 : float, optional
        Oxygen gas molar fraction. Defaults to ambient gas composition of 0.21
    P : float, optional
        Total pressure of the system. Defaults to ambient gas pressure of 1 atm/Bar.

    Returns
    -------
    tuple of length 2
        Reuturns 2 numpy arrays: Sr vacancies in mol per unit LSCF volume, delta_G in J/mol as a function of temperature given in T_range.

    Warnings
    --------
    The physical relevance of this case is shady, as vacancies resulting from SrO formation is in the bulk of the system and not near the surface. 

    """
    p_O2 = x_O2 * P
    delta_G_list = []
    V_Sr= []
    delta_E = E_LSCF_bulk_Sr_vac + E_SrO_epitax - (E_LSCF_bulk + 0.5*E_DFT_O2) + E_int
    for T in T_range:
        delta_G = E_LSCF_bulk_Sr_vac + chem_pot_SrO(T) - (E_LSCF_bulk + 0.5*(chem_pot_O2(T, E_DFT_O2, P=P))) + E_int
        if T==973: print("case 5", delta_G/ev2J_p_mol)
        delta_G_list.append(delta_G)

        K = np.exp(-delta_G/(R*T))
        #print(K)

        #DEFINITION OF THIRD DEGREE POLYNOMIAL THAT IS REACTION MECHANISM SPECIFIC
        #SHOULD BE VERIFIED EVERY TIME THE REACTION MECHANISM IS CHANGED
        N = K * np.sqrt(p_O2/P)
        a = 4+4*N
        b= 4*(x-1*N)
        c= x**2 + (1-x)*N * (1+3*x)
        d = -N * x * (1-x)**2
        solution = cubic_model(a,b,c,d)
        V_Sr.append(solution[0])
    return(np.asarray(V_Sr), np.asarray(delta_G_list))

def case6(T_range, x= 0.4):
    """Case where LSCF bulk system is used instead of slab systems. Oxygen comes from the LSCF.

    Parameters
    ----------
    T_range : numpy array
        Temperature window in Kelvins
    x : float, optional
        Initial amount of Sr present in LSCF, in mol per unit LSCF volume. Defaults to standard LSCF composition of 0.4.

    Returns
    -------
    tuple of length 2
        Returns two arrays: Sr vacacnies in moles per unit LSCF volume, delta_G in J/mol as a function of temperature given in T_range.

    Warnings
    --------
    The physical relevance of this case is shady, as vacancies resulting from SrO formation is in the bulk of the system and not near the surface. 

    """
    delta_G_list = []
    V_Sr=[]
    delta_E = E_LSCF_bulk_SrO_vac +  E_SrO_epitax - (E_LSCF_bulk) + E_int
    for T in T_range:
        delta_G = E_LSCF_bulk_SrO_vac +  chem_pot_SrO(T) - (E_LSCF_bulk) + E_int
        if T==973: print("case 6", delta_G/ev2J_p_mol)
        delta_G_list.append(delta_G)

        K = np.exp(-delta_G/(R*T))
        a = (1-(1/K))
        b = -3-x
        c = 3*x
        V_Sr.append(quadratic_model(a,b,c,x))
    return(np.asarray(V_Sr), np.asarray(delta_G_list))
