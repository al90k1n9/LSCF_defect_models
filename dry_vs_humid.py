import numpy as np
from lib.SrO_dry_air_models import case4 as dry_case
from lib.SrO_hydroxylated_models import case2 as hydrox_case
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

x = 0.4 #molar fraction of Sr
x_O2 = 0.21
P = 1 #atm
T_lower_bound = 700
T_upper_bound = 1299
T_range = np.arange(T_lower_bound,T_upper_bound,1) #K

p_O2 = x_O2 * P

fig, ax = plt.subplots(layout="constrained")
default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

ylist, delta_G_list = dry_case(T_range, x = 0.4)
ax.plot(T_range, ylist, color="red", label="dry")


ylist, delta_G, theta_list = hydrox_case(T_range)
ax.plot(T_range, ylist, color = default_colors[0], label="hydroxylated")

ax.set_xlim(T_lower_bound, T_upper_bound+1)
ax.set_ylim(0,)

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

ax.set_xlabel("T[K]")
ax.set_ylabel("$[V\'\'\'_{La}]_{eq}=x_{eq}$")

def yaxconvert(x):
    return x * 100/0.4

def yaxinvert(x):
    return x *0.4/100

ax.legend(loc="upper right", facecolor="none")


secyax = ax.secondary_yaxis("right", functions=(yaxconvert, yaxinvert))
secyax.set_ylabel("% of initial Sr content $\\frac{100 \\cdot x_{eq}}{x_0}$")
secyax.yaxis.set_minor_locator(AutoMinorLocator())

plt.show()
