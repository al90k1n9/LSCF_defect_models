import sys
sys.path.append("H:\Documents\SrO_defect_model")


from lib.chemical_potentials import *
from lib.dft_energies_0K import * 
from lib.auxilliary_functions import * 


def case1(T_range, x=0.4, p_H2O = 0.08, P = 1):
    V_Sr = np.zeros(len(T_range))
    delta_G_list = np.zeros(len(T_range))

    for index in range(len(T_range)):
        T = T_range[index]
        mu_H2O = cp_H2O(T, E_DFT_H2O=E_DFT_H2O)
        mu_SrOH2 = cp_SrOH2(T)
        delta_G = 0.5*(2*mu_SrOH2 + E_LSCF_slab_Sr_surf_O_sub_surf - 2*mu_H2O - E_LSCF_slab)
        delta_G_list[index] = delta_G
        A = 1
        B = np.exp(-delta_G/(R*T))
        C = -0.4*np.exp(-delta_G/(R*T))
        V_Sr[index] = quadratic_model(A,B,C)
    return V_Sr