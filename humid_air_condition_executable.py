import matplotlib.pyplot as plt
from lib.humid_models import *

x = 0.4 #molar fraction of Sr
x_O2 = 0.21
x_H2O = 0.08
P = 1 #atm
T_lower_bound = 700
T_upper_bound = 1400
T_range = np.arange(T_lower_bound,T_upper_bound,1) #K
#numpy imported chemical potentials, which is imported in humid models

p_O2 = x_O2 * P
p_H2O = x_H2O * P

V_Sr = case1(T_range, x, p_O2, p_H2O, P=1)
plt.plot(T_range, V_Sr)
plt.show()