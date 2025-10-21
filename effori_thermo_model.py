import numpy as np
import matplotlib.pyplot as  plt
import  os
from scipy.optimize import fsolve

local_path = os.path.dirname(os.path.abspath(__file__))
local_path += "/"

data973 = np.genfromtxt(local_path + "lib/bouwmeester_973.csv", delimiter=",", skip_header=1)

print(np.shape(data973))

#reduced model with O minus
T_ref = 700+273.15 #K
po2_ref = 0.21 #atm
#theta_o_minus_ref = 1e-4 #ref condition at T=700, po2= 0.21
gamma = 1e-5   #mol/m^2
Comax_ref = 83108 #mol/m^3
Coeq_ref = 82491 #mol/m^3
cdl = 0.2 #F/m^2

faraday = 9.648533212e4 #C/mol
R = 8.314 #J/(K.mol)

def equation4_gen(theta_o_minus, theta_o, theta_o2, K4):
    xi_eq =  gamma*faraday/cdl * theta_o_minus
    rhs = theta_o / theta_o_minus * np.exp(-faraday * xi_eq /(R*T_ref))
    return rhs-K4

def equation5_gen(theta_o_minus, theta_o, theta_o2, K5):
    theta_s =  1-(theta_o + theta_o_minus + theta_o2)
    rhs  = theta_o2 * theta_s/theta_o**2
    return rhs - K5

def equation6_gen(theta_o_minus, theta_o, theta_o2, po2, K6):
    theta_s =  1-(theta_o + theta_o_minus + theta_o2)
    rhs =  theta_s * po2 /theta_o2
    return rhs - K6

def coeq_calc(theta_o_minus, theta_o, theta_o2, K2):
    xi_eq =  gamma * faraday/cdl * theta_o_minus
    theta_s  = 1 - (theta_o_minus + theta_o + theta_o2)
    dinominator = K2*theta_s/theta_o_minus * np.exp(faraday*xi_eq/(R*T_ref))  + 1
    return Comax_ref/dinominator


po2_range= np.arange(0.1, 1.01, 0.01)
po2_range = [0.1, 0.21, 0.4, 0.8]

def elementary_model(cdl=cdl, theta_o_minus_ref = 6e-3, theta_o_ref = 1e-3, theta_o2_ref = 1e-4, gamma = 1e-5, graph=True, graph_label ="", verbose = True):
    xi_ref =  gamma * faraday/cdl * theta_o_minus_ref
    theta_s_ref = 1 - (theta_o_minus_ref + theta_o_ref + theta_o2_ref)
    #defintion of thermodynamic constants at ref condition
    K_ass = theta_o2_ref * theta_s_ref/theta_o_ref**2
    K_exc = (Comax_ref-Coeq_ref)/Coeq_ref * theta_o_minus_ref/theta_s_ref * np.exp(-faraday*xi_ref/(R*T_ref))
    K_ads = po2_ref * theta_s_ref/theta_o2_ref
    K_deion = theta_o_ref/theta_o_minus_ref * np.exp(-faraday*xi_ref/(R*T_ref))

    print(K_ass, K_exc, K_ads, K_deion)

    #Creating lists for variables that are functions of po2
    theta_o_minus_list = []
    theta_o_list = []
    theta_o2_list = []
    Coeq_list = []

    for po2 in po2_range:
        equation4 = lambda x,y, z:equation4_gen(x,y,z, K4 = K_deion)
        equation5 = lambda x, y, z:equation5_gen(x, y, z, K5 = K_ass)
        equation6 = lambda x, y, z:equation6_gen(x, y, z, po2= po2, K6 = K_ads)
        equations  = [equation4, equation5, equation6]
        initial_guess =  [1e-10, 1e-10, 1e-10] #initial guess for the surface coverages

        solution = fsolve(lambda vars:[equation(vars[0], vars[1], vars[2]) for equation in equations], initial_guess, xtol=1e-10)
        theta_o_minus_list.append(solution[0])
        theta_o_list.append(solution[1])
        theta_o2_list.append(solution[2])
        Coeq_list.append(coeq_calc(solution[0], solution[1], solution[2], K2=K_exc))
    print("po2 ", po2_range)
    print("coeq ", Coeq_list)
    print(theta_o_minus_list)
    print(theta_o_list)
    print(theta_o2_list)
    print(np.ones(len(theta_o_list))-(np.asarray(theta_o_list) +  np.asarray(theta_o_minus_list) + np.asarray(theta_o2_list)))
    fig,ax = plt.subplots(layout="constrained")
    ax.plot(po2_range,Coeq_list)


elementary_model()
plt.show()
