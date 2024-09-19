#this file needs to be completed later with the appropriate defect models.
#for now it will return the delta G values at 1000K and delta_E at 0K
#case description to be added into the github repository on the main page later

from lib.chemical_potentials import *
from lib.dft_energies_0K import * #importing all the values of dft energies
from lib.auxilliary_functions import *

def case1(T_range, x=0.4,  p_O2 = 0.21, p_H2O = 0.08, P=1):
    #see cintia mail from 19/08/2024 for case description
    #this two water molecule to form hydroxide and hydrogen gas
    #Sr from surface
    delta_E = ((E_LSCF_slab_Sr_vac_surf + 2 * E_DFT_H2 + E_SrOH2_bulk*2) - (E_LSCF_slab +  4*E_DFT_H2O))/2
    print(delta_E/ev2J_p_mol)
    V_Sr = []
    for T in T_range:
        mu_H2O = cp_H2O(T, E_DFT_H2O)
        mu_H2 = cp_H2(T, E_DFT_H2)
        p_H2 = pH2_giver(T, p_H2O, p_O2)
        delta_G = ((E_LSCF_slab_Sr_vac_surf + 2 * mu_H2 + E_SrOH2_bulk*2) - (E_LSCF_slab +  4*mu_H2O))/2

        K = np.exp(-delta_G/(R*T))
        N = K * (p_H2O/P)**2 * (P/p_H2)

        a = 4+4*N
        b= 4*(x-1*N)
        c= x**2 + (1-x)*N * (1+3*x)
        d = -N * x * (1-x)**2
        solution= cubic_model(a,b,c,d)
        V_Sr.append(solution[0])
    return V_Sr

def case2(T_range, x=0.4, p_O2 = 0.21, p_H2O = 0.08, P=1):
    #uses half oxygen gas and one water moelcule to form hydroxide
    #Sr from surface
    delta_E = ((E_LSCF_slab_Sr_vac_surf + E_SrOH2_bulk *2 ) - (2*E_DFT_H2O + E_DFT_O2 + E_LSCF_slab))/2
    print(delta_E/ev2J_p_mol)
    V_Sr = []
    for T in T_range:
        mu_H2O = cp_H2O(T, E_DFT_H2O)
        mu_O2 = cp_O2(T, E_DFT_O2)

        delta_G = ((E_LSCF_slab_Sr_vac_surf + E_SrOH2_bulk *2 ) - (2*mu_H2O + mu_O2 + E_LSCF_slab))/2
        K = np.exp(-delta_G/(R*T))
        N = K * p_H2O/P * np.sqrt(p_O2/P)

        a = 4+4*N
        b= 4*(x-1*N)
        c= x**2 + (1-x)*N * (1+3*x)
        d = -N * x * (1-x)**2
        solution= cubic_model(a,b,c,d)
        V_Sr.append(solution[0])
    return V_Sr

def case3(T_range, x=0.4, p_O2 = 0.21, p_H2O = 0.08, P=1):
    #case 1 but with both hydrogen atoms stabilised in the sr vacancy
    delta_E = ((E_LSCF_double_hydrogenated + 2* E_SrOH2_bulk) - (E_LSCF_slab + 4*E_DFT_H2O))/2
    print(delta_E/ev2J_p_mol)
    V_Sr = []
    for T in T_range:
        mu_H2O = cp_H2O(T, E_DFT_H2O)
        delta_G = ((E_LSCF_double_hydrogenated + 2* E_SrOH2_bulk) - (E_LSCF_slab + 4*mu_H2O))/2
        K = np.exp(-delta_G/(R*T))
        N = (p_H2O/P)**2 * K
        V_Sr.append(N/(1+N) * x)
    return V_Sr

def case4(T_range, x=0.4, p_O2 = 0.21, p_H2O = 0.08, P=1):
    #one hydrogen stabilised in vacancy and the other as hydrogen gas
    delta_E = ((E_LSCF_single_hydrogenated + 2* E_SrOH2_bulk + E_DFT_H2) - (E_LSCF_slab + 4*E_DFT_H2O))/2
    print(delta_E/ev2J_p_mol)
    V_Sr = []
    for T in T_range:
        mu_H2 = cp_H2(T,E_DFT_H2)
        mu_H2O = cp_H2O(T, E_DFT_H2O)
        p_H2 = pH2_giver(T, p_H2O, p_O2)
        delta_G = ((E_LSCF_single_hydrogenated + 2* E_SrOH2_bulk + mu_H2) - (E_LSCF_slab + 4*mu_H2O))/2
        K = np.exp(-delta_G/(R*T))
        N = (p_H2O/P)**2 * (P/p_H2)**(1/2) * K

        a = 1-N
        b = x+N
        c = -N*(x-x**2)
        V_Sr.append(quadratic_model(a,b,c))

    return V_Sr

def case5(T_range, x=0.4, p_O2 = 0.21, p_H2O = 0.08, P=1):
    #one water molecule and one oxygen from the material to form SrOH2
    delta_E = ((E_LSCF_slab_Sr_surf_O_sub_surf + 2*E_SrOH2_bulk) - (E_LSCF_slab + 2*E_DFT_H2O))/2
    print(delta_E/ev2J_p_mol)
    V_Sr = []
    for T in T_range:
        mu_H2O = cp_H2O(T, E_DFT_H2O)
        delta_G = ((E_LSCF_slab_Sr_surf_O_sub_surf + 2*E_SrOH2_bulk) - (E_LSCF_slab + 2*mu_H2O))/2
        K = np.exp(-delta_G/(R*T))
        N = K * p_H2O/P
        a = (1-(1/N))
        b = -3-x
        c = 3*x

        V_Sr.append(quadratic_model(a,b,c,x))

    return V_Sr