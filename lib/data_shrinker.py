import numpy as np
from lib.dft_energies_0K import ev2J_p_mol
import os

local_path = os.path.dirname(os.path.abspath(__file__))
local_path += "/"

data = np.genfromtxt(local_path + "sroh2_factsage.csv", delimiter=";")
print("input data size: ", np.shape(data))
print(data[0,0], data[-1, 0])


output_filename="SrOH2_factsage_processed.csv"
energy_threshold = 0.0001 #eV

T_list = [data[0,0]]
delta_G_list = [data[0, 2]]

for index in range (1, np.shape(data)[0]):
    if abs((data[index, 2] - delta_G_list[-1])/ev2J_p_mol)>energy_threshold:
        T_list.append(data[index,0])
        delta_G_list.append(data[index,2])


export_data = np.vstack((np.asarray(T_list), np.asarray(delta_G_list))).T
print(np.shape(export_data))
np.savetxt(local_path + output_filename, export_data, delimiter=";")
