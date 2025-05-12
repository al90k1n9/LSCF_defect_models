import matplotlib.pyplot as plt 
from lib.SrCrO4_models import *
#chemical potentials are imported in models
#numpy imported chemical potentials
from matplotlib.ticker import  AutoMinorLocator

default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

x0 = 0.4 #molar fraction of Sr
x_O2 = 0.21
x_CrO3 = 1e-10
P = 1 #atm
T_lower_bound = 700
T_upper_bound = 1200
T_range = np.linspace(T_lower_bound,T_upper_bound,1000) #K

p_O2 = x_O2 * P

fig, ax = plt.subplots(layout="constrained")
fig2, ax2 = plt.subplots(layout="constrained")
T = 800
V_Sr, delta_G, delta_G_instance,(a,b,c,d) = case1(T_range, x_CrO3=x_CrO3)
print(V_Sr[0], V_Sr[0] - V_Sr[-1])



#polynomial = []
#og_function = []
#xlist = np.linspace(0.28, 0.32, 100)
#for x in xlist:
#	polynomial.append(a*x**3 + b*x**2 + c*x + d)
#	numerator = x * (x0+2*x)**2
#	denominator = (x0-x)*(1-x0-2*x)**2 * x_CrO3 * np.sqrt(x_O2)
#	og_function.append(numerator/denominator - np.exp(-delta_G_instance/(R*T)))


#plt.plot(xlist, polynomial)
#plt.plot(xlist, og_function)
#plt.gca().axhline(y=0, color="black")
##plt.ylim(-0.1, 0.1)
#plt.show()

ax.plot(T_range, np.asarray(delta_G)/ev2J_p_mol, label="case 1")
ax2.plot(T_range, V_Sr[:,0], label="case 1 solution 1", color=default_colors[0])
ax2.plot(T_range, V_Sr[:,1], label="case 1 solution 2", color=default_colors[0], ls="dotted")
ax2.plot(T_range, V_Sr[:,2], label="case 1 solution 3", color=default_colors[0], ls = "dashed")



V_Sr, delta_G = case2(T_range, x_CrO3=x_CrO3)

ax.plot(T_range, delta_G/ev2J_p_mol, label="case 2")
ax2.plot(T_range, V_Sr, label="case 2", color=default_colors[1])


#T_range = np.genfromtxt("lib/CrO2OH2_factsage.csv", delimiter=";")[:,0]
#T_range = T_range[25:]
#print(T_range[0], T_range[-1], np.shape(T_range))
V_Sr, delta_G, mu_list = case3(T_range)

ax2.plot(T_range, V_Sr[:,0], color=default_colors[2], label="case3")
ax2.plot(T_range, V_Sr[:,1], color=default_colors[2], ls="dashed", label="case3")
ax2.plot(T_range, V_Sr[:,2], color=default_colors[2], ls="dotted", label="case3")
ax.plot(T_range, delta_G/ev2J_p_mol, color=default_colors[2], label="case 3")


ax.set_xlabel("T [K]")
ax.set_ylabel("$\Delta_rG^*(T, p)$")

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

ax.legend(facecolor = "none")

def yaxconvert(x):
    return x * 100/0.4

def yaxinvert(x):
    return x *0.4/100

secyax = ax2.secondary_yaxis("right", functions=(yaxconvert, yaxinvert))
secyax.set_ylabel("% of initial Sr content $\\frac{100 \\cdot x_{eq}}{x_0}$")
secyax.yaxis.set_minor_locator(AutoMinorLocator())

ax2.set_xlabel("T [K]")
ax2.set_ylabel("$[V_{La}''']$")

ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())

ax2.legend(facecolor="none")
ax2.set_xlim(T_lower_bound, T_upper_bound+1)
#ax2.set_ylim(0,)


fig3, ax3 = plt.subplots(layout="constrained")
ax3.plot(T_range, np.asarray(mu_list)/ev2J_p_mol)

ax3.set_xlabel("T [K]")
ax3.set_ylabel("[eV]")


ax3.set_xlim(T_lower_bound, T_upper_bound+1)


ax3.xaxis.set_minor_locator(AutoMinorLocator())
ax3.yaxis.set_minor_locator(AutoMinorLocator())

plt.show()
