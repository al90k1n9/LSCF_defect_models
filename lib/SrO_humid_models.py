import sys
sys.path.append("H:\Documents\SrO_defect_model")
from lib.chemical_potentials import *
from lib.dft_energies_0K import * #importing all the values of dft energies
from lib.auxilliary_functions import *
from H_vibration import *

def case1(T_range, x0=0.4, x_O2 = 0.21, x_H2O = 0.08, P=1):
    p_O2 = x_O2 * P
    p_H2O = x_H2O * P
    V_Sr= []
    delta_E = (E_LSCF_slab_Sr_vac_surf + 2*E_SrO_epitax + 2*E_DFT_H2- (E_LSCF_slab + 2*E_DFT_H2O))/2 + E_int
    print(delta_E/ev2J_p_mol, " delta_E of case 1")
    delta_G_range = []
    p_H2_list = []
    for T in T_range:
        p_H2 = pH2_giver(T, p_H2O, p_O2)
        p_H2_list.append(p_H2)

        delta_G = (E_LSCF_slab_Sr_vac_surf + 2*chem_pot_SrO(T) + 2*chem_pot_H2(T, E_DFT_H2, P=P)- (E_LSCF_slab + 2*chem_pot_H2O(T, E_DFT_H2O, P=P)))/2 + E_int
        delta_G_range.append(delta_G)
        if T == 1000: print("delta G at T = 1000", delta_G/ev2J_p_mol)
        K = np.exp(-delta_G/(R*T))

        #DEFINITION OF THIRD DEGREE POLYNOMIAL THAT IS REACTION MECHANISM SPECIFIC
        #SHOULD BE VERIFIED EVERY TIME THE REACTION MECHANISM IS CHANGED
        N = K * p_H2O/p_H2 #notice that the total pressure cancels out in the case
        a = 4+4*N
        b= 4*(x0-1*N)
        c= x0**2 + (1-x0)*N * (1+3*x0)
        b = -4*N
        c = (1-x0)*N * (1+3*x0)
        d = -N * x0 * (1-x0)**2
        solution= cubic_model(a,b,c,d)
        #print(x0_minus- x0_plus)
        #print(equation(x0_minus), equation(x0_plus))
        V_Sr.append(solution[0])
        #print(solution[0])
    return (np.asarray(V_Sr), np.asarray(delta_G_range), np.asarray(p_H2_list))

def case2(T_range, x0=0.4, x_H2O = 0.08, P=1):
    p_H2O = x_H2O * P
    V_Sr= []
    delta_E = (2*E_SrO_epitax + E_LSCF_double_hydrogenated - (E_LSCF_slab+ 2 * E_DFT_H2O)) / 2 + E_int
    print(delta_E/ev2J_p_mol, " delta_E of case 2")
    delta_G_range = []
    for T in T_range:
        delta_G = (2*chem_pot_SrO(T) + E_LSCF_double_hydrogenated - (E_LSCF_slab+ 2 *chem_pot_H2O(T, E_DFT_H2O, P=P))) / 2  + E_int + 2*OH_bond_vibration
        delta_G_range.append(delta_G)
        if T == 1000: print("case 2: delta G at T = 1000", delta_G/ev2J_p_mol)
        K = np.exp(-delta_G/(R*T)) 
        N = K * p_H2O/P 
        V_Sr.append(N/(1+N)*x0)

    return (np.asarray(V_Sr), np.asarray(delta_G_range))

def case3(T_range, x0=0.4, x_O2 = 0.21, x_H2O = 0.08, P=1):
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
    return (V_Sr, delta_G_range)

def case4(T_range, x0=0.4, x_O2 = 0.21, x_H2O = 0.08, P=1):
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
    return V_Sr
