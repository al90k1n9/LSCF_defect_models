#this file needs to be completed later with the appropriate defect models.
#for now it will return the delta G values at 1000K and delta_E at 0K
#case description to be added into the github repository on the main page later

from lib.chemical_potentials import *
from lib.dft_energies_0K import * #importing all the values of dft energies
from lib.auxilliary_functions import *

def case1(T, p_O2 = 0.21, p_H2O = 0.08, P=1):
    #see cintia mail from 19/08/2024 for case description
    #this two water molecule to form hydroxide and hydrogen gas
    #Sr from surface
    delta_E = ((E_LSCF_slab_Sr_vac_surf + 2 * E_DFT_H2 + E_SrOH2_bulk*2) - (E_LSCF_slab +  4*E_DFT_H2O))/2

    mu_H2O = cp_H2O(T, E_DFT_H2O)
    mu_H2 = cp_H2(T, E_DFT_H2)

    delta_G = ((E_LSCF_slab_Sr_vac_surf + 2 * mu_H2 + E_SrOH2_bulk*2) - (E_LSCF_slab +  4*mu_H2O))/2
    return (delta_E/ev2J_p_mol, delta_G/ev2J_p_mol)

def case2(T, p_O2 = 0.21, p_H2O = 0.08, P=1):
    #uses half oxygen gas and one water moelcule to form hydroxide
    #Sr from surface
    delta_E = ((E_LSCF_slab_Sr_vac_surf + E_SrOH2_bulk *2 ) - (2*E_DFT_H2O + E_DFT_O2 + E_LSCF_slab))/2
    
    mu_H2O = cp_H2O(T, E_DFT_H2O)
    mu_O2 = cp_O2(T, E_DFT_O2)

    delta_G = ((E_LSCF_slab_Sr_vac_surf + E_SrOH2_bulk *2 ) - (2*mu_H2O + mu_O2 + E_LSCF_slab))/2
    return (delta_E/ev2J_p_mol, delta_G/ev2J_p_mol)

def case3(T, p_O2 = 0.21, p_H2O = 0.08, P=1):
    #case 1 but with hydrogen stabilised in the sr vacancy
    delta_E = ((E_LSCF_double_hydrogenated + 2* E_SrOH2_bulk) - (E_LSCF_slab + 4*E_DFT_H2O))/2

    mu_H2O = cp_H2O(T, E_DFT_H2O)
    delta_G = ((E_LSCF_double_hydrogenated + 2* E_SrOH2_bulk) - (E_LSCF_slab + 4*mu_H2O))/2
    return (delta_E/ev2J_p_mol, delta_G/ev2J_p_mol)

def case4(T, p_O2 = 0.21, p_H2O = 0.08, P=1):
    delta_E = ((E_LSCF_single_hydrogenated + 2* E_SrOH2_bulk + E_DFT_H2) - (E_LSCF_slab + 4*E_DFT_H2O))/2
    mu_H2 = cp_H2(T,E_DFT_H2)
    mu_H2O = cp_H2O(T, E_DFT_H2O)
    delta_G = ((E_LSCF_single_hydrogenated + 2* E_SrOH2_bulk + mu_H2) - (E_LSCF_slab + 4*mu_H2O))/2
    return (delta_E/ev2J_p_mol, delta_G/ev2J_p_mol)

def case5(T, p_O2 = 0.21, p_H2O = 0.08, P=1):
    delta_E = ((E_LSCF_slab_Sr_surf_O_sub_surf + 2*E_SrOH2_bulk) - (E_LSCF_slab + 2*E_DFT_H2O))/2
    mu_H2O = cp_H2O(T, E_DFT_H2O)
    delta_G = ((E_LSCF_slab_Sr_surf_O_sub_surf + 2*E_SrOH2_bulk) - (E_LSCF_slab + 2*mu_H2O))/2
    return (delta_E/ev2J_p_mol, delta_G/ev2J_p_mol)