import numpy as np
from lib.SrO_dry_air_models import case4 as dry_case4
from lib.SrO_dry_air_models import case1 as dry_case1
from lib.SrO_hydroxylated_models import case2 as hydrox_case
from lib.SrO_hydroxylated_models import case1 as hydrox_case1
from lib.SrO_humid_models import case1 as humid_case1
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

x = 0.4 #molar fraction of Sr
x_O2 = 0.21
P = 1 #atm
T_lower_bound = 900
T_upper_bound = 1100
T_range = np.arange(T_lower_bound,T_upper_bound,1) #K

p_O2 = x_O2 * P

fig, ax = plt.subplots(layout="constrained")
default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

ylist, delta_G_list = dry_case4(T_range, x = 0.4)
ax.plot(T_range, ylist, label="dry case 4")

ylist, delta_G_list = dry_case1(T_range, x = 0.4)
ax.plot(T_range, ylist, label="dry case 1")

ylist, delta_G, theta_list = hydrox_case(T_range)
ax.plot(T_range, ylist,  label="hydroxylated 2H stabilisied")

ylist, delta_G, theta_list = hydrox_case1(T_range)
ax.plot(T_range, ylist,  label="hydroxylated H2")

ylist, delta_G, theta_list = humid_case1(T_range)
ax.plot(T_range, ylist,  label="humid H2")

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
