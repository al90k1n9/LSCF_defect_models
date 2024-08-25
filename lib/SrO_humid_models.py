from lib.chemical_potentials import *
from lib.dft_energies_0K import * #importing all the values of dft energies
from lib.auxilliary_functions import *
from H_vibration import vibrational_correction_term

def case1(T_range, x=0.4, p_O2 = 0.21, p_H2O = 0.08, P=1):
    V_Sr= []
    delta_E = (E_LSCF_slab_Sr_vac_surf + 2*E_SrO + 2*E_DFT_H2- (E_LSCF_slab + 2*E_DFT_H2O))/2
    print(delta_E/ev2J_p_mol, " delta_E of case 1")
    for T in T_range:
        mu_H2O = cp_H2O(T, E_DFT_H2O)
        p_H2 = pH2_giver(T, p_H2O, p_O2)
        mu_H2 = cp_H2(T, E_DFT_H2)

        delta_G = (E_LSCF_slab_Sr_vac_surf + 2*E_SrO + 2*mu_H2- (E_LSCF_slab + 2*mu_H2O))/2
        if T == 1000: print("delta G at T = 1000", delta_G/ev2J_p_mol)
        K = np.exp(-delta_G/(R*T))

        #DEFINITION OF THIRD DEGREE POLYNOMIAL THAT IS REACTION MECHANISM SPECIFIC
        #SHOULD BE VERIFIED EVERY TIME THE REACTION MECHANISM IS CHANGED
        N = K * p_H2O/p_H2 #notice that the total pressure cancels out in the case
        a = 4+4*N
        b= 4*(x-1*N)
        c= x**2 + (1-x)*N * (1+3*x)
        d = -N * x * (1-x)**2
        solution= cubic_model(a,b,c,d)
        #print(x0_minus- x0_plus)
        #print(equation(x0_minus), equation(x0_plus))
        V_Sr.append(solution[0])
        #print(solution[0])
    return V_Sr

def case2(T_range, x=0.4, p_H2O = 0.08, P=1):
    V_Sr= []
    delta_E = (2*E_SrO + E_LSCF_double_hydrogenated - (E_LSCF_slab+ 2 * E_DFT_H2O)) / 2
    for T in T_range:
        mu_H2O = cp_H2O(T, E_DFT_H2O)
        delta_G = (2*E_SrO + E_LSCF_double_hydrogenated - (E_LSCF_slab+ 2 * mu_H2O)) / 2 - vibrational_correction_term(T)
        K = np.exp(-delta_G/(R*T)) 
        N = K * p_H2O/P 
        V_Sr.append(N/(1+N)*x)

    return V_Sr

def case3(T_range, x=0.4, p_O2 = 0.21, p_H2O = 0.08, P=1):
    V_Sr= []
    delta_E = (E_LSCF_single_hydrogenated + 2*E_SrO + E_DFT_H2- (E_LSCF_slab + 2*E_DFT_H2O))/2
    print(delta_E/ev2J_p_mol, " delta_E of case 3")
    for T in T_range:
        mu_H2O = cp_H2O(T, E_DFT_H2O)
        p_H2 = pH2_giver(T, p_H2O, p_O2)
        mu_H2 = cp_H2(T, E_DFT_H2)

        delta_G = (E_LSCF_single_hydrogenated + 2*E_SrO + mu_H2- (E_LSCF_slab + 2*mu_H2O))/2
        #if T == 1000: print("delta G at T = 1000", delta_G/ev2J_p_mol)
        K = np.exp(-delta_G/(R*T))

        #DEFINITION OF THIRD DEGREE POLYNOMIAL THAT IS REACTION MECHANISM SPECIFIC
        #SHOULD BE VERIFIED EVERY TIME THE REACTION MECHANISM IS CHANGED
        N = K * p_H2O/P * np.sqrt(P/p_H2)
        a = 1-N
        b = x+N
        c = -N*(x-x**2)
        V_Sr.append(quadratic_model(a,b,c))
    return V_Sr

def case4(T_range, x=0.4, p_O2 = 0.21, p_H2O = 0.08, P=1):
    V_Sr= []
    delta_E = (E_LSCF_slab_Sr_vac_bulk + E_DFT_H2 + E_SrO) - (E_LSCF_slab + E_DFT_H2O)
    print(delta_E/ev2J_p_mol, " delta_E of case 4")
    for T in T_range:
        mu_H2O = cp_H2O(T, E_DFT_H2O)
        p_H2 = pH2_giver(T, p_H2O, p_O2)
        mu_H2 = cp_H2(T, E_DFT_H2)

        delta_G = (E_LSCF_slab_Sr_vac_bulk + mu_H2 + E_SrO) - (E_LSCF_slab + mu_H2O)
        if T == 1000: print("delta G at T = 1000", delta_G/ev2J_p_mol)
        K = np.exp(-delta_G/(R*T))

        #DEFINITION OF THIRD DEGREE POLYNOMIAL THAT IS REACTION MECHANISM SPECIFIC
        #SHOULD BE VERIFIED EVERY TIME THE REACTION MECHANISM IS CHANGED
        N = K * p_H2O/p_H2 #notice that the total pressure cancels out in the case
        a = 4+4*N
        b= 4*(x-1*N)
        c= x**2 + (1-x)*N * (1+3*x)
        d = -N * x * (1-x)**2
        solution= cubic_model(a,b,c,d)
        #print(x0_minus- x0_plus)
        #print(equation(x0_minus), equation(x0_plus))
        V_Sr.append(solution[0])
        #print(solution[0])
    return V_Sr
