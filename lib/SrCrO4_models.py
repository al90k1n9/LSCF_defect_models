from lib.chemical_potentials import *
from lib.dft_energies_0K import *
from lib.auxilliary_functions import *

def case1(T_range, x0=0.4, x_O2=0.21, x_CrO3 = 1e-3, P=1):
	p_O2 = P*x_O2
	V_Sr_list=[]
	delta_G_list = []
	delta_E = (E_LSCF_slab_Sr_vac_surf + 2*E_DFT_SrCrO4 - 2*E_DFT_CrO3 - E_DFT_O2 - E_LSCF_slab)/2
	for T in T_range:
		mu_CrO3 = cp_CrO3(T, E_DFT_CrO3, P=1)
		mu_O2 = cp_O2(T, E_DFT_O2, P=1)
		delta_G = (E_LSCF_slab_Sr_vac_surf + 2*E_DFT_SrCrO4 - 2*mu_CrO3 - mu_O2 - E_LSCF_slab)/2
		delta_G_list.append(delta_G)
		N = x_CrO3 * np.sqrt(x_O2)  * np.exp(-delta_G/(R*T))

		a = -4*N-4
		b = 4*N - 4*x0
		c = - (N * (1-x0)**2 + 4*x0* (1-x0)* N + x0**2) 
		d = N*x0*(1-x0)**2
		solution = cubic_model(a,b,c,d)

		V_Sr_list.append(np.sort(solution))
	return (np.asarray(V_Sr_list), delta_G_list, delta_G, (a,b,c,d))

def case2(T_range, x0=0.4, x_CrO3 = 1e-3, P = 1):
	V_Sr_list=[]
	delta_G_list = []
	delta_E = (E_LSCF_slab_Sr_surf_O_sub_surf + 2*E_DFT_SrCrO4 - E_DFT_O2 - E_LSCF_slab)/2
	for T in T_range:
		mu_CrO3 = cp_CrO3(T, E_DFT_CrO3, P=1)
		delta_G = (E_LSCF_slab_Sr_surf_O_sub_surf + 2*E_DFT_SrCrO4 - 2*mu_CrO3 - E_LSCF_slab)/2
		delta_G_list.append(delta_G)
		N = x_CrO3 * np.exp(-delta_G/(R*T))

		a = 1 - 1/N
		b = -(3+x0)
		c = 3*x0 

		solution = quadratic_model(a,b,c)

		V_Sr_list.append(solution)
	return (np.asarray(V_Sr_list), np.asarray(delta_G_list))

def case3(T_range, x0= 0.4, x_chr_hyd = 1e-3, x_O2 = 0.21, x_H2O = 0.08, P=1):
    V_Sr_list= [] 
    delta_G_list = []
    delta_E = (E_LSCF_slab_Sr_vac_surf + 2*E_DFT_SrCrO4 + 2*E_DFT_H2 - (E_LSCF_slab + 2* E_DFT_CrO2OH2))*1/2
    mu_list = []
    for T in T_range:
        mu_CrO2OH2 = cp_CrO2OH2(T, P=1)
        mu_H2 = cp_H2(T, E_DFT_H2, P=1)
        mu_list.append(mu_CrO2OH2-E_DFT_CrO2OH2)
        p_H2O = P * x_H2O
        p_O2 = x_O2 * P
        p_H2 = pH2_giver(T, p_H2O, p_O2)
        delta_G = (E_LSCF_slab_Sr_vac_surf + 2*E_DFT_SrCrO4 + 2*mu_H2 - (E_LSCF_slab + 2* mu_CrO2OH2))*1/2
        delta_G_list.append(delta_G)
        N = x_chr_hyd / (p_H2/P) * np.exp(-delta_G/(R*T))
        a = -4*N-4
        b = 4*N - 4*x0
        c = - (N * (1-x0)**2 + 4*x0* (1-x0)* N + x0**2) 
        d = N*x0*(1-x0)**2
        solution = cubic_model(a,b,c,d)
        V_Sr_list.append(solution)
    return (np.asarray(V_Sr_list), np.asarray(delta_G_list), mu_list)
