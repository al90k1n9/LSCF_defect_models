import numpy as np


N_avagadro = 6.0223*10**23
ev2J = 1.60219*10**(-19)
ev2J_p_mol = ev2J*N_avagadro
Ha2eV = 27.2114
OH_bond_vibration = 0.351/2 * ev2J_p_mol #J/mol

acell_SrO_slab = 5.16756734063135E-10 #meters #this is with standard SrO lattice parameter
acell_LSCF_slab = 7.74799694590908E-10 #meters
acell_LSCF_slab_expaned =  7.84E-10 #meters

#atoms in isolation
E_Sr =  -8.46031896932468E+02 #eV
E_O = -4.32296299818140E+02 #eV

E_Sr *= ev2J_p_mol
E_O *= ev2J_p_mol

E_SrO = -1.28749128492132E+03 #eV
E_SrO = E_SrO * ev2J_p_mol #in J/mol

E_SrO_epitax = -1.28729253918189E+03 #eV
E_SrO_epitax *= ev2J_p_mol

E_SrO_expanded_substrate = -1.28719164856155E+03 #eV
E_SrO_expanded_substrate *= ev2J_p_mol


E_SrO_slab_expanded_substrate = -1.28711321076406E+04 #eV
E_SrO_slab_expanded_substrate *= ev2J_p_mol



E_DFT_O2 = -874.671000 #eV
E_DFT_O2 = E_DFT_O2 * ev2J_p_mol #J/mol

E_DFT_H2O = -4.71505500161390E+02 #eV
E_DFT_H2O *= ev2J_p_mol #J/mol

E_DFT_SrCrO4 = -1.99412155436326E+04/4
E_DFT_SrCrO4 *= ev2J_p_mol

E_DFT_CrO3 = -3.69256819208978E+03
E_DFT_CrO3 *= ev2J_p_mol


E_DFT_CrO2OH2 = -4.16229857029462E+03
E_DFT_CrO2OH2 *= ev2J_p_mol

E_DFT_H2 = -3.17817399567036E+01  #in ev without zero point energy corrention at T=0K
E_DFT_H2 = E_DFT_H2 * ev2J_p_mol #in J/mol

E_SrOH2_bulk = -7.04103776E+03 * ev2J_p_mol/4 #the system has four molecules and therefore the factor 1/4

E_LSCF_slab_Sr_vac_surf = -9.37015381643914E+04 #in eV
E_LSCF_slab_Sr_vac_surf = E_LSCF_slab_Sr_vac_surf * ev2J_p_mol #in J/mol

E_LSCF_slab = -9.54034366543809E+04 #eV
E_LSCF_slab = E_LSCF_slab * ev2J_p_mol #in J/mol

E_LSCF_SrO_interface = -1.211567739277E+05 #eV
#the system contains 20 SrO unit cells and the LSCF slab
E_LSCF_SrO_interface *= ev2J_p_mol #J/mol

E_SrO_slab = -1.2873269133E+04 #eV SrO slab containing 10 unit SrO cells
E_SrO_slab *= ev2J_p_mol #J/mol

E_SrO_slab_epitax = -1.287173093E+04 #eV lattice parameter has been elaraged to match that of LSCF
E_SrO_slab_epitax *= ev2J_p_mol 

gamma_SrO = (E_SrO_slab - 10* E_SrO)/(2*acell_SrO_slab**2) #J/mol/m**2
gamma_SrO_epitax = (E_SrO_slab_epitax - 10* E_SrO_epitax)/(2 * 0.5* acell_LSCF_slab**2) #J/mol/m**2
gamma_SrO_expanded_substrate = (E_SrO_slab_expanded_substrate - 10* E_SrO_expanded_substrate)/(2 * 0.5* acell_LSCF_slab**2) #J/mol/m**2



adhesion_work = (E_LSCF_slab + 20 * E_SrO_epitax + 2* gamma_SrO_epitax * acell_LSCF_slab**2 - E_LSCF_SrO_interface)/(2*acell_LSCF_slab**2) #J/mol/m**2


E_int = (acell_LSCF_slab/2)**2/2 * (2*gamma_SrO_epitax - adhesion_work) #J/mol
E_int_expanded_substrate = (acell_LSCF_slab/2)**2/2 * (2*gamma_SrO_expanded_substrate - adhesion_work) #J/mol


E_LSCF_slab_Sr_vac_bulk = -94551.83065 #eV
E_LSCF_slab_Sr_vac_bulk = E_LSCF_slab_Sr_vac_bulk * ev2J_p_mol #J/mol

E_LSCF_slab_Sr_surf_O_sub_surf = -9.28256221726180E+04 #eV
E_LSCF_slab_Sr_surf_O_sub_surf = E_LSCF_slab_Sr_surf_O_sub_surf * ev2J_p_mol #J/mol

E_LSCF_slab_SrO_bulk = -9.41146616226594E+04 #eV
E_LSCF_slab_SrO_bulk = E_LSCF_slab_SrO_bulk * ev2J_p_mol #J/mol

E_LSCF_slab_SrO_bulk = -9.41146616226594E+04 #eV
E_LSCF_slab_SrO_bulk = E_LSCF_slab_SrO_bulk * ev2J_p_mol #J/mol

E_LSCF_bulk_Sr_vac = -8.93195943147177E+04 # eV
E_LSCF_bulk_Sr_vac = E_LSCF_bulk_Sr_vac * ev2J_p_mol #J/mol

E_LSCF_bulk_SrO_vac = -8.88822924941026E+04 #eV
E_LSCF_bulk_SrO_vac = E_LSCF_bulk_SrO_vac * ev2J_p_mol #J/mol

E_LSCF_bulk = -9.01712962373595E+04 #eV
E_LSCF_bulk = E_LSCF_bulk * ev2J_p_mol #J/mol

E_LSCF_hydroxilated = -9.63495122688081E+04 #eV
E_LSCF_hydroxilated *= ev2J_p_mol



E_LSCF_bulk_hydrogenated = -9.45693723052703E+04 #eV
E_LSCF_bulk_hydrogenated *= ev2J_p_mol #J/mol

E_LSCF_double_hydrogenated = -9.37728973049342E+04
E_LSCF_double_hydrogenated *= ev2J_p_mol

E_LSCF_single_hydrogenated = -9.37371667638564E+04
E_LSCF_single_hydrogenated *= ev2J_p_mol 



#=========================================================================================================================
#Vibrational properties of hydrogen bonds in single hydrogenated LSCF slab
#=========================================================================================================================
single_hydrogenation_configs = np.asarray([-3444.77767328802, -3444.7358570430, -3444.7357648852, -3444.7565289008, -3444.7728438132, -3444.7638908796, -3444.7737535218]) * Ha2eV * ev2J_p_mol
#======================================================================================================================================================================================================================================



zpe_H2O = 0.56 * ev2J_p_mol
zpe_H2 = 0.27 * ev2J_p_mol
zpe_O2 = 0.098 * ev2J_p_mol
zpe_CrO3 = 0.25 * ev2J_p_mol

E_ads = (E_LSCF_hydroxilated - (E_LSCF_slab + 2* E_DFT_H2O))/2  #J/mol

double_hydrogenation_energy = (E_LSCF_double_hydrogenated - (E_LSCF_slab_Sr_vac_surf + 2 * E_DFT_H2))/2

single_hydrogenation_energy = (E_LSCF_single_hydrogenated-(E_LSCF_slab_Sr_vac_surf + E_DFT_H2))/2

second_hydrogenation_energy = (E_LSCF_double_hydrogenated - (E_LSCF_single_hydrogenated + E_DFT_H2))/2


if __name__ == "__main__":
    print(single_hydrogenation_energy/ev2J_p_mol, " single hydrogenation_energy")
    print(double_hydrogenation_energy/ev2J_p_mol, " double hydrogenation_energy")
    print(second_hydrogenation_energy/ev2J_p_mol, " second hydrogenation_energy")
    print((second_hydrogenation_energy+single_hydrogenation_energy)/ev2J_p_mol, " sum of first and second")
    
    print("interface energy contribution in eV")
    print(E_int/ev2J_p_mol, " E_int")
    print((E_SrO - (E_Sr + E_O))/ev2J_p_mol, " cohesive energy SrO")
    print((E_SrO_epitax + E_int - (E_Sr + E_O))/ev2J_p_mol, " energy to put SrO on lscf surface")
    print((E_SrO_expanded_substrate + E_int - (E_Sr + E_O))/ev2J_p_mol, "energy to put SrO on lscf substrate, with lscf exapnded due to high temperatures")
    print((E_SrO-E_SrO_epitax)/ev2J_p_mol, " difference in energy hen lattice contracted")
    print(E_int_expanded_substrate/ev2J_p_mol)
    
    
    print("SrO free surface energies in eV/[1x1]")
    print(gamma_SrO/ev2J_p_mol*acell_SrO_slab**2/2, " SrO lattice parameter at room temperature")
    print(gamma_SrO_epitax/ev2J_p_mol*acell_LSCF_slab**2*0.25, " SrO lattice parameter matching LSCF lattice parameter from calculations")
    print(gamma_SrO_expanded_substrate/ev2J_p_mol*acell_LSCF_slab_expaned**2*0.25, " LSCF lattice parameter expanded for temperature effects and then SrO lattice parameter matched\n")
    
    print((E_SrO_expanded_substrate-E_SrO_epitax)/ev2J_p_mol, " difference in unit cell energy")
