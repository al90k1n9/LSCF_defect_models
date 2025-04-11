import sys
sys.path.append("H:\Documents\SrO_defect_model")
from lib.chemical_potentials import *
from lib.dft_energies_0K import * #importing all the values of dft energies
from lib.auxilliary_functions import *
from H_vibration import *

sro_vibration_data = np.genfromtxt("./lib/vibrational_correction_sro.csv",delimiter=" ")
sro_vibration_data[:,1] *= ev2J_p_mol * 0#to convert everything in J/mol units
T_data = sro_vibration_data[:,0]


def case1(T_range, x=0.4, p_O2 = 0.21, p_H2O = 0.08, P=1):
    V_Sr= []
    delta_E = (E_LSCF_slab_Sr_vac_surf + 2*E_SrO_epitax + 2*E_DFT_H2- (E_LSCF_slab + 2*E_DFT_H2O))/2 + E_int
    print(delta_E/ev2J_p_mol, " delta_E of case 1")
    delta_G_range = []
    for T in T_range:
        mu_H2O = cp_H2O(T, E_DFT_H2O, P=P)
        p_H2 = pH2_giver(T, p_H2O, p_O2)
        mu_H2 = cp_H2(T, E_DFT_H2, P=P)

        T_index_vib_data = np.where(T_data == T)
        delta_G = (E_LSCF_slab_Sr_vac_surf + 2*(E_SrO_epitax + float(sro_vibration_data[T_index_vib_data, 1])) + 2*(mu_H2 + zpe_H2)- (E_LSCF_slab + 2*(mu_H2O + zpe_H2O)))/2 + E_int
        delta_G_range.append(delta_G)
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
    return (V_Sr, delta_G_range)

def case2(T_range, x=0.4, p_H2O = 0.08, P=1):
    V_Sr= []
    delta_E = (2*E_SrO_epitax + E_LSCF_double_hydrogenated - (E_LSCF_slab+ 2 * E_DFT_H2O)) / 2 + E_int
    print(delta_E/ev2J_p_mol, " delta_E of case 2")
    delta_G_range = []
    for T in T_range:
        mu_H2O = cp_H2O(T, E_DFT_H2O, P=P)
        T_index_vib_data = np.where(T_data == T)
        delta_G = (2*(E_SrO_epitax + float(sro_vibration_data[T_index_vib_data, 1])) + E_LSCF_double_hydrogenated - (E_LSCF_slab+ 2 *(mu_H2O + zpe_H2O))) / 2 + vibrational_correction_term(T) + E_int
        if len(delta_G_range)>=1:    
            if delta_G * delta_G_range[-1]<0: print(T, " inversion temperature in K")
        delta_G_range.append(delta_G)
        if T == 1000: print("case 2: delta G at T = 1000", delta_G/ev2J_p_mol)
        K = np.exp(-delta_G/(R*T)) 
        N = K * p_H2O/P 
        V_Sr.append(N/(1+N)*x)

    return (V_Sr, delta_G_range)

def case3(T_range, x=0.4, p_O2 = 0.21, p_H2O = 0.08, P=1):
    V_Sr= []
    delta_E = (E_LSCF_single_hydrogenated + 2*E_SrO_epitax + E_DFT_H2- (E_LSCF_slab + 2*E_DFT_H2O))/2 + E_int
    print(delta_E/ev2J_p_mol, " delta_E of case 3")
    delta_G_range = []
    for T in T_range:
        mu_H2O = cp_H2O(T, E_DFT_H2O, P=P)
        p_H2 = pH2_giver(T, p_H2O, p_O2)
        mu_H2 = cp_H2(T, E_DFT_H2, P=P)

        T_index_vib_data = np.where(T_data == T)
        delta_G = (E_LSCF_single_hydrogenated + 2*(E_SrO_epitax + float(sro_vibration_data[T_index_vib_data, 1])) + mu_H2 + zpe_H2- (E_LSCF_slab + 2*(mu_H2O + zpe_H2O)))/2 + vibrational_correction_term(T)/2 + E_int
        delta_G_range.append(delta_G)
        if T == 1000: print("case 3: delta G at T = 1000", delta_G/ev2J_p_mol)
        K = np.exp(-delta_G/(R*T))

        #DEFINITION OF THIRD DEGREE POLYNOMIAL THAT IS REACTION MECHANISM SPECIFIC
        #SHOULD BE VERIFIED EVERY TIME THE REACTION MECHANISM IS CHANGED
        N = K * p_H2O/P * np.sqrt(P/p_H2)
        a = 1-N
        b = x+N
        c = -N*(x-x**2)
        V_Sr.append(quadratic_model(a,b,c))
    return (V_Sr, delta_G_range)

def case4(T_range, x=0.4, p_O2 = 0.21, p_H2O = 0.08, P=1):
    V_Sr= []
    delta_E = (E_LSCF_slab_Sr_vac_bulk + E_DFT_H2 + E_SrO_epitax) - (E_LSCF_slab + E_DFT_H2O) + E_int
    print(delta_E/ev2J_p_mol, " delta_E of case 4")
    delta_G_range = []
    for T in T_range:
        mu_H2O = cp_H2O(T, E_DFT_H2O, P=P)
        p_H2 = pH2_giver(T, p_H2O, p_O2)
        mu_H2 = cp_H2(T, E_DFT_H2, P=P)

        T_index_vib_data = np.where(T_data == T)
        delta_G = (E_LSCF_slab_Sr_vac_bulk + mu_H2 + zpe_H2 + (E_SrO_epitax + float(sro_vibration_data[T_index_vib_data, 1]))) - (E_LSCF_slab + mu_H2O + zpe_H2O) + E_int
        delta_G_range.append(delta_G)
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
    return (V_Sr, delta_G_range)
