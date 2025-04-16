from lib.chemical_potentials import *
from lib.dft_energies_0K import *
from lib.auxilliary_functions import *

def case1(T_range, x0=0.4, x_O2=0.21, P=1, x_CrO3 = 1e-3):
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

		V_Sr_list.append(solution[0])
	return (V_Sr_list, delta_G_list, delta_G, (a,b,c,d))
