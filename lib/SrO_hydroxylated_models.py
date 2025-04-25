from lib.chemical_potentials import *
from lib.dft_energies_0K import * #importing all the values of dft energies
from lib.auxilliary_functions import *



sro_vibration_data = np.genfromtxt("./lib/vibrational_correction_sro.csv",delimiter=" ")
sro_vibration_data[:,1] +=0
sro_vibration_data[:,1] *= 0* ev2J_p_mol #to convert everything in J/mol units
T_data = sro_vibration_data[:,0]

def case1(T_range, x=0.4, p_O2 = 0.21, p_H2O = 0.08, P=1):
    V_Sr= []
    delta_E = (E_LSCF_slab_Sr_vac_surf + 2*E_SrO_epitax + 2*E_DFT_H2- E_LSCF_hydroxilated )/2 + E_int
    print(delta_E/ev2J_p_mol, " of case 1")
    delta_G_list = []
    theta_list = []
    for T in T_range:
        p_H2 = pH2_giver(T, p_H2O, p_O2)
        mu_H2 = cp_H2(T, E_DFT_H2, P=P)

        T_index_vib_data = np.where(T_data == T)
        delta_G = (E_LSCF_slab_Sr_vac_surf + 2*(E_SrO_epitax + float(sro_vibration_data[T_index_vib_data, 1])) + 2*(mu_H2+zpe_H2)- E_LSCF_hydroxilated )/2 + E_int
        delta_G_list.append(delta_G)
        theta = surface_coverage_H2O(T,p_H2O/P, E_ads, P)
        N = theta/(1-theta) * np.exp(-delta_G/(R*T)) * P/p_H2

        a = 4+4*N
        b= 4*(x-1*N)
        c= x**2 + (1-x)*N * (1+3*x)
        d = -N * x * (1-x)**2
        solution= cubic_model(a,b,c,d)
        V_Sr.append(solution[0])
        theta_list.append(theta)
    return (V_Sr, np.asarray(delta_G_list), theta_list)

def case2(T_range, x=0.4, p_H2O = 0.08, P=1):
    V_Sr= []
    delta_E = (E_LSCF_double_hydrogenated + 2*E_SrO_epitax - E_LSCF_hydroxilated )/2 + E_int
    print(delta_E/ev2J_p_mol, " of case 2")
    
    delta_G_list = []
    theta_list = []
    for T in T_range:
        T_index_vib_data = np.where(T_data == T)
        delta_G = delta_E + float(sro_vibration_data[T_index_vib_data, 1]) 
        delta_G_list.append(delta_G)
        theta = surface_coverage_H2O(T,p_H2O/P, E_ads, P)
        N = theta/(1-theta) * np.exp(-delta_G/(R*T))
        V_Sr.append(N/(1+N)*x)
        theta_list.append(theta)
    return (np.asarray(V_Sr), np.asarray(delta_G_list), theta_list)

def case3(T_range, x=0.4, p_O2 = 0.21, p_H2O = 0.08, P=1):
    V_Sr= []
    delta_E = (E_LSCF_single_hydrogenated + 2*E_SrO_epitax + E_DFT_H2 - E_LSCF_hydroxilated )/2 + E_int
    print(delta_E/ev2J_p_mol, " of case 3")
    delta_G_list = []
    theta_list= []
    for T in T_range:
        p_H2 = pH2_giver(T, p_H2O, p_O2)
        mu_H2 = cp_H2(T, E_DFT_H2, P=P)

        T_index_vib_data = np.where(T_data == T)
        delta_G = (E_LSCF_single_hydrogenated + 2*(E_SrO_epitax + float(sro_vibration_data[T_index_vib_data, 1])) + mu_H2 + zpe_H2 - E_LSCF_hydroxilated )/2 + E_int
        delta_G_list.append(delta_G)

        theta = surface_coverage_H2O(T,p_H2O/P, E_ads, P)
        N = theta/(1-theta) * np.exp(-delta_G/(R*T)) * np.sqrt(P/p_H2)

        a = 1-N
        b = x+N
        c = -N*(x-x**2)
        V_Sr.append(quadratic_model(a,b,c))
        theta_list.append(theta)
    return (V_Sr, np.asarray(delta_G_list), theta_list)
