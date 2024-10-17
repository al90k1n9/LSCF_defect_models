import numpy as np


N_avagadro = 6.0223*10**23
ev2J = 1.60219*10**(-19)
ev2J_p_mol = ev2J*N_avagadro
Ha2eV = 27.2114

E_SrO = -1.28749128492132E+03 #eV
E_SrO = E_SrO * ev2J_p_mol #in J/mol

E_DFT_O2 = -874.671000 #eV
E_DFT_O2 = E_DFT_O2 * ev2J_p_mol #J/mol

E_DFT_H2O = -4.71505500161390E+02 #eV
E_DFT_H2O *= ev2J_p_mol #J/mol

E_DFT_H2 = -3.17817399567036E+01  #in ev without zero point energy corrention at T=0K
E_DFT_H2 = E_DFT_H2 * ev2J_p_mol #in J/mol

E_SrOH2_bulk = -7.04103776E+03 * ev2J_p_mol/4 #the system has four molecules and therefore the factor 1/4

E_LSCF_slab_Sr_vac_surf = -9.37015381643914E+04 #in eV
E_LSCF_slab_Sr_vac_surf = E_LSCF_slab_Sr_vac_surf * ev2J_p_mol #in J/mol

E_LSCF_slab = -9.54034366543809E+04
E_LSCF_slab = E_LSCF_slab * ev2J_p_mol #in J/mol

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

E_ads = (E_LSCF_hydroxilated - (E_LSCF_slab + 2* E_DFT_H2O))/2 #J/mol

E_LSCF_bulk_hydrogenated = -9.45693723052703E+04 #eV
E_LSCF_bulk_hydrogenated *= ev2J_p_mol #J/mol

E_LSCF_double_hydrogenated = -9.37728973049342E+04
E_LSCF_double_hydrogenated *= ev2J_p_mol

E_LSCF_single_hydrogenated = -9.37371667638564E+04 * ev2J_p_mol 


double_hydrogenation_energy = (E_LSCF_double_hydrogenated - (E_LSCF_slab_Sr_vac_surf + 2 * E_DFT_H2))/2

singe_hydrogenation_energy = (E_LSCF_single_hydrogenated-(E_LSCF_slab_Sr_vac_surf + E_DFT_H2))/2

#=========================================================================================================================
#Vibrational properties of hydrogen bonds in single hydrogenated LSCF slab
#=========================================================================================================================
single_hydrogenation_configs = np.asarray([-3444.77767328802, -3444.7358570430, -3444.7357648852, -3444.7565289008, -3444.7728438132, -3444.7638908796, -3444.7737535218]) * Ha2eV * ev2J_p_mol
#======================================================================================================================================================================================================================================
