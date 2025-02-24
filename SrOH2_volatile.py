from lib.SrOH2_volatile_models import *
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator



x = 0.4 #molar fraction of Sr
x_O2 = 0.21
x_H2O = 0.08
P = 1 #atm
T_lower_bound = 700 #K
T_upper_bound = 1400 #K


data = np.genfromtxt("./lib/sroh2_factsage.csv", delimiter=";")
T_range = data[:,0]


V_Sr = case1(T_range)


fig, ax = plt.subplots(layout="constrained")
ax.plot(T_range, V_Sr)

ax.set_xlabel("T [K]")
ax.set_ylabel("V_Sr")

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

plt.show()
