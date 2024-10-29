from lib.chemical_potentials import *
from lib.dft_energies_0K import * #importing all the values of dft energies
from lib.auxilliary_functions import *

def case1(T_range, x=0.4, p_O2 = 0.21, p_H2O = 0.08, P=1):
    V_Sr= []
    delta_E = (E_LSCF_slab_Sr_vac_surf + 2*E_SrO + 2*E_DFT_H2- E_LSCF_hydroxilated )/2
    print(delta_E/ev2J_p_mol, " of case 1")
    for T in T_range:
        p_H2 = pH2_giver(T, p_H2O, p_O2)
        mu_H2 = cp_H2(T, E_DFT_H2, P=P)

        delta_G = (E_LSCF_slab_Sr_vac_surf + 2*E_SrO + 2*mu_H2- E_LSCF_hydroxilated )/2
        theta = surface_coverage_H2O(T,p_H2O/P, E_ads, P)
        N = theta/(1-theta) * np.exp(-delta_G/(R*T)) * P/p_H2

        a = 4+4*N
        b= 4*(x-1*N)
        c= x**2 + (1-x)*N * (1+3*x)
        d = -N * x * (1-x)**2
        solution= cubic_model(a,b,c,d)
        V_Sr.append(solution[0])
    return V_Sr

def case2(T_range, x=0.4, p_H2O = 0.08, P=1):
    V_Sr= []
    delta_E = (E_LSCF_double_hydrogenated + 2*E_SrO - E_LSCF_hydroxilated )/2
    print(delta_E/ev2J_p_mol, " of case 2")
    delta_G = delta_E
    for T in T_range:
        theta = surface_coverage_H2O(T,p_H2O/P, E_ads, P)
        N = theta/(1-theta) * np.exp(-delta_G/(R*T))
        V_Sr.append(N/(1+N)*x)
    return V_Sr

def case3(T_range, x=0.4, p_O2 = 0.21, p_H2O = 0.08, P=1):
    V_Sr= []
    delta_E = (E_LSCF_single_hydrogenated + 2*E_SrO + E_DFT_H2- E_LSCF_hydroxilated )/2
    print(delta_E/ev2J_p_mol, " of case 3")
    for T in T_range:
        p_H2 = pH2_giver(T, p_H2O, p_O2)
        mu_H2 = cp_H2(T, E_DFT_H2, P=P)

        delta_G = (E_LSCF_single_hydrogenated + 2*E_SrO + mu_H2- E_LSCF_hydroxilated )/2

        theta = surface_coverage_H2O(T,p_H2O/P, E_ads, P)
        N = theta/(1-theta) * np.exp(-delta_G/(R*T)) * np.sqrt(P/p_H2)

        a = 1-N
        b = x+N
        c = -N*(x-x**2)
        V_Sr.append(quadratic_model(a,b,c))
    return V_Sr
