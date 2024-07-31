
#delta mu is defined as follows => mu(T,p) = E_DFT(T=0, p~0) + delta_mu(T,p)
#it's not the same delta mu as done in the derivation
import numpy as np

N_avagadro = 6.0223*10**23
ev2J = 1.60219*10**(-19)
ev2J_p_mol = ev2J*N_avagadro
R = 8.314 #J.K.mol-1

def cp_O2(T, E_DFT_O2, p=1):
    #chemical potential of a pure O2 system. There is no notion of partial pressure here.
    T_0 = 298 #K
    cp = 4.83*10**-23 #K J/atom
    P_0 = 1 #bar
    delta_h0 = 8700 #J.mol-1
    s0 = 205 #J.mol-1.K-1

    delta_mu_O2 = delta_h0 + cp* N_avagadro*(T-T_0) - T*s0 - T*cp*N_avagadro*np.log(T/T_0) + T*R * np.log(p/P_0) #J/mol

    mu_O2 = E_DFT_O2 + delta_mu_O2
    return mu_O2

def cp_H2O(T, E_DFT_H2O, P=1):
    T_0 = 298 #K
    cp = 33.6 / N_avagadro #K J/atom
    P_0 = 1 #bar
    delta_h0 = 9905 #J.mol-1
    s0 = 188.835 #J.mol-1.K-1

    delta_mu_H2O = delta_h0 + cp* N_avagadro*(T-T_0) - T*s0 - T*cp*N_avagadro*np.log(T/T_0) + T*R * np.log(P/P_0) #J/mol

    mu_H2O = E_DFT_H2O + delta_mu_H2O
    return mu_H2O


def cp_H2(T, E_DFT_H2, P=1):
    T_0 = 298 #K
    cp = 28.8 / N_avagadro #K J/atom
    P_0 = 1 #bar
    delta_h0 = 8468 #J.mol-1
    s0 = 130.68 #J.mol-1.K-1
    #delta mu is defined as follows => mu(T,p) = E_DFT(T=0, p~0) + delta_mu(T,p)
    #it's not the same delta mu as done in the derivation
    delta_mu_H2 = delta_h0 + cp* N_avagadro*(T-T_0) - T*s0 - T*cp*N_avagadro*np.log(T/T_0) + T*R * np.log(P/P_0) #J/mol
    #print(delta_mu_H2, "delta_mu_H2 in J/mol")
    mu_H2 = E_DFT_H2 + delta_mu_H2
    return mu_H2


