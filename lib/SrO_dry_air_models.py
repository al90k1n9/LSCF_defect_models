from lib.chemical_potentials import *
from lib.dft_energies_0K import * #importing all the values of dft energies
from lib.auxilliary_functions import *



def case1(T_range, x=0.4, p_O2 = 0.21, P=1):
    V_Sr= []
    delta_G_list = []
    delta_E = (E_LSCF_slab_Sr_vac_surf + 2 * E_SrO - (E_LSCF_slab + E_DFT_O2))/2
    for T in T_range:
        mu_O2 = cp_O2(T, E_DFT_O2, P=P)
        delta_G = (E_LSCF_slab_Sr_vac_surf + 2 * E_SrO - (E_LSCF_slab + mu_O2))/2
        if T==1000: print("case 1", delta_G/ev2J_p_mol)
        K = np.exp(-delta_G/(R*T))
        #print(K)

        #DEFINITION OF THIRD DEGREE POLYNOMIAL THAT IS REACTION MECHANISM SPECIFIC
        #SHOULD BE VERIFIED EVERY TIME THE REACTION MECHANISM IS CHANGED
        N = K * np.sqrt(p_O2/P)
        a = 4+4*N
        b= 4*(x-1*N)
        c= x**2 + (1-x)*N * (1+3*x)
        d = -N * x * (1-x)**2
        solution= cubic_model(a,b,c,d)
        #print(x0_minus- x0_plus)
        #print(equation(x0_minus), equation(x0_plus))
        V_Sr.append(solution[0])
        delta_G_list.append(delta_G)
        #print(solution[0])
    return (V_Sr,delta_G_list)

def case2(T_range, x=0.4, p_O2 = 0.21, P=1):
    delta_G_list =[]
    V_Sr= []
    delta_E = E_LSCF_slab_Sr_vac_bulk + E_SrO - (E_LSCF_slab + 0.5*E_DFT_O2)
    for T in T_range:
        mu_O2 = cp_O2(T, E_DFT_O2, P=P)

        delta_G = E_LSCF_slab_Sr_vac_bulk + E_SrO - (E_LSCF_slab + 0.5*mu_O2)
        if T==1000: print("case 2", delta_G/ev2J_p_mol)
        delta_G_list.append(delta_G)
        K = np.exp(-delta_G/(R*T))
        #print(K)

        #DEFINITION OF THIRD DEGREE POLYNOMIAL THAT IS REACTION MECHANISM SPECIFIC
        #SHOULD BE VERIFIED EVERY TIME THE REACTION MECHANISM IS CHANGED
        N = K * np.sqrt(p_O2/P)
        a = 4+4*N
        b= 4*(x-1*N)
        c= x**2 + (1-x)*N * (1+3*x)
        d = -N * x * (1-x)**2
        solution= cubic_model(a,b,c,d)
        V_Sr.append(solution[0])
    return (V_Sr,delta_G_list)

def case3(T_range, x=0.4):
    delta_G_list = []
    V_Sr=[]
    delta_E = E_LSCF_slab_Sr_surf_O_sub_surf /2 +  E_SrO - (E_LSCF_slab/2)
    for T in T_range:
        delta_G = delta_E
        if T==1000: print("case 3", delta_G/ev2J_p_mol)
        delta_G_list.append(delta_G)

        K = np.exp(-delta_G/(R*T))
        a = (1-(1/K))
        b = -3-x
        c = 3*x

        V_Sr.append(quadratic_model(a,b,c,x))
    return (V_Sr,delta_G_list)

def case4(T_range, x=0.4):
    delta_G_list = []
    V_Sr=[]
    delta_E = E_LSCF_slab_SrO_bulk +  E_SrO - (E_LSCF_slab)
    for T in T_range:
        delta_G = delta_E
        if T==963: print("case 4", delta_G/ev2J_p_mol)
        delta_G_list.append(delta_G)

        K = np.exp(-delta_G/(R*T))
        a = (1-(1/K))
        b = -3-x
        c = 3*x

        V_Sr.append(quadratic_model(a,b,c,x))
    return (V_Sr,delta_G_list)

def case5(T_range, x=0.4, p_O2 = 0.21, P=1):
    delta_G_list = []
    V_Sr= []
    delta_E = E_LSCF_bulk_Sr_vac + E_SrO - (E_LSCF_bulk + 0.5*E_DFT_O2)
    for T in T_range:
        mu_O2 = cp_O2(T, E_DFT_O2, P=P)
        delta_G = E_LSCF_bulk_Sr_vac + E_SrO - (E_LSCF_bulk + 0.5*mu_O2)
        if T==1000: print("case 5", delta_G/ev2J_p_mol)
        delta_G_list.append(delta_G)

        K = np.exp(-delta_G/(R*T))
        #print(K)

        #DEFINITION OF THIRD DEGREE POLYNOMIAL THAT IS REACTION MECHANISM SPECIFIC
        #SHOULD BE VERIFIED EVERY TIME THE REACTION MECHANISM IS CHANGED
        N = K * np.sqrt(p_O2/P)
        a = 4+4*N
        b= 4*(x-1*N)
        c= x**2 + (1-x)*N * (1+3*x)
        d = -N * x * (1-x)**2
        solution = cubic_model(a,b,c,d)
        V_Sr.append(solution[0])
    return (V_Sr,delta_G_list)

def case6(T_range, x= 0.4):
    delta_G_list = []
    V_Sr=[]
    delta_E = E_LSCF_bulk_SrO_vac +  E_SrO - (E_LSCF_bulk)
    for T in T_range:
        delta_G = delta_E
        if T==1000: print("case 6", delta_G/ev2J_p_mol)
        delta_G_list.append(delta_G)

        K = np.exp(-delta_G/(R*T))
        a = (1-(1/K))
        b = -3-x
        c = 3*x
        V_Sr.append(quadratic_model(a,b,c,x))
    return (V_Sr,delta_G_list)
