N_avagadro = 6.0223*10**23
ev2J = 1.60219*10**(-19)
ev2J_p_mol = ev2J*N_avagadro

E_LSCF_slab_Sr_vac_surf = -9.37015381643914E+04 #in eV
E_LSCF_slab_Sr_vac_surf = E_LSCF_slab_Sr_vac_surf * ev2J_p_mol #in J/mol

E_SrO = -1.28749128492132E+03 #eV
E_SrO = E_SrO * ev2J_p_mol #in J/mol

E_LSCF_slab = -9.54034366543809E+04
E_LSCF_slab = E_LSCF_slab * ev2J_p_mol #in J/mol

E_DFT_O2 = -874.671000 #eV
E_DFT_O2 = E_DFT_O2 * ev2J_p_mol #J/mol

E_DFT_H2O = -4.71505500161390E+02 #eV
E_DFT_H2O *= ev2J_p_mol #J/mol

E_DFT_H2 = -3.17817399567036E+01  #in ev without zero point energy corrention at T=0K
E_DFT_H2 = E_DFT_H2 * ev2J_p_mol #in J/mol

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
