import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy.optimize import curve_fit

N_avagadro = 6.0223*10**23
ev2J = 1.60219*10**(-19)
ev2J_p_mol = ev2J*N_avagadro

E = -1.28749128669295E+03 #eV"
E95 = -1.28743754684416E+03 #eV
E9 = -1.28727928605548E+03 #eV
E105 = -1.28748574120615E+03 #eV
E11 =  -1.28745019842608E+03
E1025 = -1.28749388576600E+03
E0925 = -1.28740253525102E+03
E0975 = -1.28747636466546E+03
E099 = -1.28748690220602E+03


data = np.genfromtxt("./lib/SrO_vibrational_free_energy/thermo.dat")
data_0p9V = np.genfromtxt("./lib/SrO_vibrational_free_energy/thermo_dat_0.9V.dat")
data_1p1V = np.genfromtxt("./lib/SrO_vibrational_free_energy/thermo_dat_1.1V.dat")
data_1p05V = np.genfromtxt("./lib/SrO_vibrational_free_energy/thermo_dat_1.05V.dat")
data_0p95V = np.genfromtxt("./lib/SrO_vibrational_free_energy/thermo_dat_0.95V.dat")
data_0p925V = np.genfromtxt("./lib/SrO_vibrational_free_energy/thermo_dat_0.925V.dat")
data_0p975V = np.genfromtxt("./lib/SrO_vibrational_free_energy/thermo_dat_0.975V.dat")
data_1p025V = np.genfromtxt("./lib/SrO_vibrational_free_energy/thermo_dat_1.025V.dat")
data_0p99V = np.genfromtxt("./lib/SrO_vibrational_free_energy/thermo_dat_0.99V.dat")
#T, F, E,S, C, omega
print(np.shape(data))
print(np.shape(data_0p9V))
print(np.shape(data_1p05V))
print(np.shape(data_0p95V))
print(np.shape(data_1p1V))
print(np.shape(data_0p925V))
print(np.shape(data_0p975V))
print(np.shape(data_1p025V))

data[:,1:4] *= 0.5
data_0p9V[:,1:4] *= 0.5
data_1p1V[:,1:4] *= 0.5
data_0p95V[:,1:4] *= 0.5
data_1p05V[:,1:4] *= 0.5
data_1p025V[:,1:4] *= 0.5
data_0p925V[:,1:4] *= 0.5
data_0p975V[:,1:4] *= 0.5


fig3, ax3 = plt.subplots(layout="constrained")
volume_list = np.asarray([0.88171, 0.925, 0.943562, 0.975, 1, 1.025, 1.05013, 1.09316])
energy_vib = np.asarray([data_0p9V[0,2], data_0p925V[0,2], data_0p95V[0,2], data_0p975V[0,2], data[0,2], data_1p025V[0,2], data_1p05V[0,2], data_1p1V[0,2]])
energy_0K =  np.asarray([E9, E0925, E95, E0975, E, E1025, E105, E11])


ax3.plot(volume_list, energy_0K-E, marker="o")

ax3.set_xlabel("V [#V$_{0}$]")
ax3.set_ylabel("static lattice energy [eV]")
ax3.xaxis.set_minor_locator(AutoMinorLocator())
ax3.yaxis.set_minor_locator(AutoMinorLocator())

print(np.min(energy_0K-E))


fig4, ax4 = plt.subplots(layout="constrained", figsize=(8,6))


T_range = np.arange(1,1300,150)
T_len = np.shape(T_range)[0]
V_min_index = 2
V_max_index = -1
EpF = np.zeros((T_len, 8))
eq_V = np.zeros(T_len)
F_eq = np.zeros(T_len)

for index in range(T_len):
    T = T_range[index]
    EpF[index, 0] = E9+data_0p9V[T, 1]/ev2J_p_mol
    EpF[index, 1] = E0925+data_0p925V[T,1]/ev2J_p_mol
    EpF[index, 2] = E95 + data_0p95V[T, 1]/ev2J_p_mol
    EpF[index, 3] = E0975 + data_0p975V[T, 1]/ev2J_p_mol
    EpF[index, 4] = E+data[T,1]/ev2J_p_mol
    EpF[index, 5] = E1025+data_1p025V[T,1]/ev2J_p_mol
    EpF[index, 6] = E105+data_1p05V[T, 1]/ev2J_p_mol
    EpF[index, 7] = E11 + data_1p1V[T, 1]/ev2J_p_mol
    #print(np.asarray([E9+data_0p9V[T, 1], E95 + data_0p95V[T, 1], E+data[T,1], E105+data_1p05V[T, 1], E11 + data_1p1V[T, 1]]))

def quadratic_fit(x, a,b,c):
    return a*x**2 + b* x + c

xlist = volume_list[V_min_index:V_max_index]
art_xlist = np.arange(volume_list[V_min_index], volume_list[V_max_index-1], 0.001)
for index in range(T_len):
    ylist = EpF[index,V_min_index:V_max_index]-E
    popt, pcov = curve_fit(quadratic_fit, xlist, ylist, [0,0,0])  
    ax4.plot(xlist , ylist, marker="o", markerfacecolor="None", label="T="+str(T_range[index])+"K")
    ax4.plot(art_xlist, quadratic_fit(art_xlist, popt[0], popt[1], popt[2]), color="black", ls="dashed")
    xmin = -popt[1]/(2*popt[0])
    F_min = quadratic_fit(xmin, popt[0], popt[1], popt[2])
    eq_V[index] = xmin
    F_eq[index] = F_min
    ax4.plot(xmin, F_min, marker = 'o', markersize=5, color="black")
ax4.legend(loc="upper center", facecolor="none", ncol=5, bbox_to_anchor=(0.5, 1.15))

ax4.set_ylabel("F$_{vib}$ + E$_{static}(V)$ - E$_{static}$(V=V$_0$) [eV]")
ax4.set_xlabel("Volume [#V$_0$]")


fig, ax = plt.subplots(layout="constrained")
ax.set_xlabel("T [K]")
ax.set_ylabel("cell parameter [${\AA}$]")

ax.set_xlim(0, T_range[-1]+1)

ax.plot(T_range, (eq_V * 3.65**3)**(1/3), markerfacecolor="None", color= "black")
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())


fig2, ax2 = plt.subplots(layout="constrained")
ax2.set_xlabel("T [K]")
ax2.set_ylabel("F$_{vib}$ + E$_{static}$(V) - E$_{static}$(V=V0) [eV]")
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())

ax2.plot(T_range, F_eq, label="F$_{vib}$(Veq)")
#ax2.plot(T_range, data[:-1, 1]/ev2J_p_mol, label="F$_{vib}$ at V0")

ax2.legend(loc = "upper right", facecolor="none")

print(np.shape(T_range), np.shape(F_eq))
np.savetxt("vibrational_correction_sro.csv", np.column_stack((T_range, F_eq)), delimiter=" ")



fig5, ax5 = plt.subplots(layout="constrained")
ax5.set_xlabel("T[K]")
ax5.set_ylabel("difference in vibrtional free energy [eV]")

#ax5.plot(T_range, F_eq-data[:-1,1]/ev2J_p_mol)

fig6, ax6 = plt.subplots(layout="constrained")
ax6.plot(T_range, F_eq)
ax6.set_xlabel("T[K]")
ax6.set_ylabel("F$_{vib}$ at V$_{eq}$ [eV]")
ax6.xaxis.set_minor_locator(AutoMinorLocator())
ax6.yaxis.set_minor_locator(AutoMinorLocator())

ax6.set_xlim(T_range[0]-1, T_range[-1]+1)



#fig.savefig("./figs/sro_phonon_lattice_expansion.svg", format="svg", dpi=300, transparent=True)
#fig3.savefig("./figs/sro_phonon_static_energy.svg", format="svg", dpi=300, transparent=True)
#fig6.savefig("./figs/sro_phonon_free_energy.svg", format="svg", dpi=300, transparent=True)
fig4.savefig("./figs/sro_phonon_feq_veq.svg", format="svg", dpi= 300, transparent=True)
plt.show()
