from lib.chemical_potentials import *
from lib.dft_energies_0K import * #importing all the values of dft energies
from lib.auxilliary_functions import *

def case1(T_range, x=0.4, p_O2 = 0.21, p_H2O = 0.08, P=1):
    #case 1
    V_Sr= []
    delta_H = (E_LSCF_slab_Sr_vac_surf + 2*E_SrO + 2*E_DFT_H2- (E_LSCF_slab + 2*E_DFT_H2O))/2
    print("reaction enthalpy: ", delta_H/ev2J_p_mol)
    for T in T_range:
        mu_H2O = cp_H2O(T, E_DFT_O2)
        p_H2 = pH2_giver(T, p_H2O, p_O2)
        mu_H2 = cp_H2(T, E_DFT_H2)

        #case 1
        delta_G = (E_LSCF_slab_Sr_vac_surf + 2*E_SrO + 2*mu_H2- (E_LSCF_slab + 2*E_DFT_H2O))/2
        K = np.exp(-delta_G/(R*T))
        #print(K)

        #DEFINITION OF THIRD DEGREE POLYNOMIAL THAT IS REACTION MECHANISM SPECIFIC
        #SHOULD BE VERIFIED EVERY TIME THE REACTION MECHANISM IS CHANGED
        N = K * np.sqrt(p_H2O/p_H2) #notice that the total pressure cancels out in the case
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