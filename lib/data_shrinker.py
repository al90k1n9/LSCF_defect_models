import numpy as np
from dft_energies_0K import ev2J_p_mol

data = np.genfromtxt("lib/CrO2OH2_factsage.csv", delimiter=";")

output_filename="CrO2OH2_factsage_processed.csv"
energy_threshold = 0.0001 #eV

T_list = [data[0,0]]
delta_G_list = [data[0, 2]]

for index in range (1, np.shape(data)[0]):
    if (data[index, 2] - delta_G_list[-1])/ev2J_p_mol>energy_threshold:
        T_list.append(data[index,0])
        delta_G_list.append(data[index,2])


export_data = np.vstack((np.asarray(T_list), np.asarray(delta_G_list))).T
print(np.shape(export_data))
np.savetxt("lib"+output_filename, export_data, delimiter=";")
