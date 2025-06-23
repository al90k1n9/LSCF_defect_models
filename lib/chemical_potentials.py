"""Chemical potentials used in the defect modelling. It is important to keep in that zpe should be added without failure.

"""

#delta mu is defined as follows => mu(T,p) = E_DFT(T=0, p~0) + delta_mu(T,p)
#it's not the same delta mu as done in the derivation
import sys
import numpy as np
from lib.dft_energies_0K import E_DFT_H2O, E_DFT_CrO3, E_SrO, E_SrO_epitax, zpe_H2O, zpe_CrO3 ,zpe_O2, zpe_H2
import os

local_path = os.path.dirname(os.path.abspath(__file__))
local_path += "/"

N_avagadro = 6.0223*10**23
ev2J = 1.60219*10**(-19)
ev2J_p_mol = ev2J*N_avagadro
R = 8.31446261815 #J.K.mol-1

data_SrOH2 = np.genfromtxt(local_path + "SrOH2_factsage_processed.csv", delimiter=";") #SrO (solid) +H2O (gas) gives SrOH2 (gas)
data_CrO2OH2 = np.genfromtxt(local_path + "CrO2OH2_factsage_processed.csv", delimiter=";") #SrO (solid) +H2O (gas) gives SrOH2 (gas)
data_SrO_vibration = np.genfromtxt(local_path + "sro_phonon_vibrational_free_energy.csv", delimiter =";") #results from dfpt

def linear_interpolator(x, x_data, y_data):
    """Function to interpolate data.

    Parameters
    ----------
    x : float
        x for which you want to determine y
    x_data : numpy array
        data
    y_data : numpy array
        data

    Returns
    -------
    float
        interpolated value
    
    Warnings
    --------
    This function cannot be moved to auxilliary functions due to circular imports. auxilliary modules imports this module to get some for surface_coverage. Maybe put surface coverage in hydroxylated_case and move this to auxilliary function for a neater structure.

    """
    assert (x>=x_data[0] and x<=x_data[-1]), "requested data not in interpolation domain" + str(x)
    if x in x_data:
        x_index = np.where(x_data == x)
        return float(y_data[x_index])
    else:
        for index in range(0, len(x_data)):
            if x_data[index] < x:
                slope = (y_data[index] - y_data[index-1])/(x_data[index] - x_data[index-1])
                intercept = y_data[index] - slope*x_data[index]
                value = slope*x+intercept
                return float(value)

def chem_pot_O2(T, E_DFT_O2, P=1):
    """Chemical potential of oxygen gas at reference condition.

    Parameters
    ----------
    T : float or numpy array
        Temperature or temperature window in Kelvins
    E_DFT_O2 : float
        Oxygen molecule energy in J/mol from DFT.
    P : float, optional
        Total pressure of the system. Defaults to ambient total pressure of 1 atm.

    Returns
    -------
    same type as T
        chemical potential of oxygen gas in J/mol with zpe included.

    """
    #chemical potential of a pure O2 system. There is no notion of partial pressure here.
    T_0 = 298 #K
    cp = 4.83*10**-23 #K J/atom
    P_0 = 1 #bar
    delta_h0 = 8700 #J.mol-1
    s0 = 205 #J.mol-1.K-1

    delta_mu_O2 = delta_h0 + cp* N_avagadro*(T-T_0) - T*s0 - T*cp*N_avagadro*np.log(T/T_0) + T*R * np.log(P/P_0) #J/mol

    mu_O2 = E_DFT_O2 + delta_mu_O2 + zpe_O2
    return mu_O2

def chem_pot_H2O(T, E_DFT_H2O, P=1):
    """Chemical potential of water vapour at reference condition.

    Parameters
    ----------
    T : float or numpy array
        Temperature or temperature window in Kelvins
    E_DFT_O2 : float
        Water molecule energy in J/mol from DFT.
    P : float, optional
        Total pressure of the system. Defaults to ambient total pressure of 1 atm.

    Returns
    -------
    same type as T
        chemical potential of water vapour in J/mol with zpe included.

    """
    T_0 = 298 #K
    cp = 33.6 / N_avagadro #K J/atom
    P_0 = 1 #bar
    delta_h0 = 9905 #J.mol-1
    s0 = 188.835 #J.mol-1.K-1

    delta_mu_H2O = delta_h0 + cp* N_avagadro*(T-T_0) - T*s0 - T*cp*N_avagadro*np.log(T/T_0) + T*R * np.log(P/P_0) #J/mol

    mu_H2O = E_DFT_H2O + delta_mu_H2O + zpe_H2O
    return mu_H2O


def chem_pot_H2(T, E_DFT_H2, P=1):
    """Chemical potential of hydrogen gas at reference condition.

    Parameters
    ----------
    T : float or numpy array
        Temperature or temperature window in Kelvins
    E_DFT_O2 : float
        Hydrogen molecule energy in J/mol from DFT.
    P : float, optional
        Total pressure of the system. Defaults to ambient total pressure of 1 atm.

    Returns
    -------
    same type as T
        chemical potential of hydrogen gas in J/mol with zpe included.

    """
    T_0 = 298 #K
    cp = 28.8 / N_avagadro #K J/atom
    P_0 = 1 #bar
    delta_h0 = 8468 #J.mol-1
    s0 = 130.68 #J.mol-1.K-1
    #delta mu is defined as follows => mu(T,p) = E_DFT(T=0, p~0) + delta_mu(T,p)
    #it's not the same delta mu as done in the derivation
    delta_mu_H2 = delta_h0 + cp* N_avagadro*(T-T_0) - T*s0 - T*cp*N_avagadro*np.log(T/T_0) + T*R * np.log(P/P_0) #J/mol
    #print(delta_mu_H2, "delta_mu_H2 in J/mol")
    mu_H2 = E_DFT_H2 + delta_mu_H2 + zpe_H2
    return mu_H2


def chem_pot_SrOH2(T, P=1, data = None):
    """Determinaton of SrOH2 chemical_potential from formation Gibbs free energies taken from free version of Factsage

    Parameters
    ----------
    T : float
        Temperature in Kelvin
    P : float, optional
        Total pressure. Defaults to 1 ambient gas pressure
    data : None
        Initialised to formation Gibbs free energy from factsage. Defaulted to none for sphinx.

    Returns
    -------
    float
        Chemical potential of SrOH2 in J/mol

    """
    data = data_SrOH2
    delta_G_sroh2 = linear_interpolator(T,data[:,0], data[:,1])
    mu_H2O = chem_pot_H2O(T, E_DFT_H2O=E_DFT_H2O)
    mu_SrOH2 = delta_G_sroh2 + mu_H2O + E_SrO + R*T*np.log(P) #J/mol
    return mu_SrOH2



def chem_pot_CrO3(T, E_DFT_CrO3, P=1):
    """Chemical potential of CrO3 gas.

    Parameters
    ----------
    T : numpy array or float
        Temperature or temperature window in Kelvins
    E_DFT_CrO3 : float
        CrO3 molecule energy in J/mol from DFT.
    P : float, optional
        Total pressure. Defaults to ambient total pressure of 1 atm/Bar

    Returns
    -------
    same type as T
        Chemical potential of CrO3 gas.

    Warnings
    --------
    ZPE of this molecule needs to be calculated.

    """
    #TODO zpe
    T_0 = 298 #K
    cp = 56.025/N_avagadro #K J/atom
    P_0 = 1 #bar
    delta_h0 = 0 #J.mol-1
    s0 = 266.178 #J.mol-1.K-1

    delta_mu_CrO3 = delta_h0 + cp* N_avagadro*(T-T_0) - T*s0 - T*cp*N_avagadro*np.log(T/T_0) + T*R * np.log(P/P_0) #J/mol

    mu_CrO3 = E_DFT_CrO3 + delta_mu_CrO3
    return mu_CrO3


def chem_pot_CrO2OH2(T, P=1, data = None):
    """Chemical potential of CrO2OH2 from formation Gibbs energies taken from free version of Factsage.

    Parameters
    ----------
    T : float
        Temperature in Kelvin
    P : float, optional
        Total pressure that defaults to ambient total pressure of 1 atm/Bar.
    data : None
        Reinitialised to formation Gibbs free energies of CrO2OH2. Defaults to None for sphinx.

    Returns
    -------
    float
        Chemical potential of CrO2OH2

    Warnings
    --------
    ZPE of CrO3 will affect this chemical potential too.

    """
    #TODO zpe of CrO3 affects this as well
    data = data_CrO2OH2
    delta_Gf_CrO2OH2 = linear_interpolator(T, data[:,0],data[:,1])
    mu_CrO2OH2 = delta_Gf_CrO2OH2 + chem_pot_H2O(T, E_DFT_H2O=E_DFT_H2O) + chem_pot_CrO3(T, E_DFT_CrO3=E_DFT_CrO3) + R*T*np.log(P) #J/mol
    return mu_CrO2OH2


def chem_pot_SrO(T, E_DFT_SrO=E_SrO_epitax, data = None, vib=0):
    """Chemical potential of SrO

    Parameters
    ----------
    T : float
        Temperature in Kelvin
    E_DFT_SrO : float, optional
        Static lattice energy of SrO from DFT in J/mol
    data : None
        Reinitialised with vibrational Helmholtz free energy in eV 
    vib : bool, optional
        Bool to include vibrational properties or not. Defaults to 0, to not include.

    Returns
    -------
    float
        Chemical potential of SrO

    """
    data = data_SrO_vibration
    if vib:
        F_vib = linear_interpolator(T, data[:,0], data[:,1]) * ev2J_p_mol
        return F_vib + E_DFT_SrO
    return E_DFT_SrO
