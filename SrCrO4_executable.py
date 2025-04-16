import matplotlib.pyplot as plt 
from lib.SrCrO4_models import *
#chemical potentials are imported in models
#numpy imported chemical potentials
from matplotlib.ticker import  AutoMinorLocator

x0 = 0.4 #molar fraction of Sr
x_O2 = 0.21
x_CrO3 = 1e-3
P = 1 #atm
T_lower_bound = 700
T_upper_bound = 1299
T_range = np.arange(T_lower_bound,T_upper_bound,1) #K

p_O2 = x_O2 * P

#fig, ax = plt.subplots(layout="constrained")
#fig2, ax2 = plt.subplots(layout="constrained")
T = 800
V_Sr, delta_G, delta_G_instance,(a,b,c,d) = case1([T])

polynomial = []
og_function = []
xlist = np.linspace(0, 0.4, 100)
for x in xlist:
	polynomial.append(a*x**3 + b*x**2 + c*x + d)
	numerator = x * (x0+2*x)**2
	denominator = (x0-x)*(1-x0-2*x)**2 * x_CrO3 * np.sqrt(x_O2)
	og_function.append(numerator/denominator - np.exp(-delta_G_instance/(R*T)))


plt.plot(xlist, polynomial)
plt.plot(xlist, og_function)
plt.gca().axhline(y=0, color="black")
plt.ylim(-0.1, 0.1)
plt.show()
#ax.plot(T_range, np.asarray(delta_G)/ev2J_p_mol, label="case 1")
#ax2.plot(T_range, V_Sr, label="case 1")
#
#
#ax.set_xlabel("T [K]")
#ax.set_ylabel("$\Delta_rG^*(T, p)$")
#
#ax.xaxis.set_minor_locator(AutoMinorLocator())
#ax.yaxis.set_minor_locator(AutoMinorLocator())
#
#ax.legend(facecolor = "none")
#
#ax2.set_xlabel("T [K]")
#ax2.set_ylabel("$[V_{La}''']$")
#
#ax2.xaxis.set_minor_locator(AutoMinorLocator())
#ax2.yaxis.set_minor_locator(AutoMinorLocator())
#
#ax2.legend(facecolor="none")


plt.show()
