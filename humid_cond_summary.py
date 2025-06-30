from lib.SrO_humid_models import case2 as humid
from lib.SrO_humid_models import case1 as humid_case1
from lib.SrO_hydroxylated_models import case2 as hydrox
from lib.SrOH2_volatile_models import case1 as sroh2_case
from matplotlib.ticker import FormatStrFormatter as fsf

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os

local_path = os.path.dirname(os.path.abspath(__file__))
local_path += "/"

x = 0.4 #molar fraction of Sr
x_O2 = 0.21
P = 1 #atm
T_lower_bound = 700
T_upper_bound = 1100
T_range = np.arange(T_lower_bound,T_upper_bound,1) #K


fig, ax = plt.subplots(layout="constrained")
default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
axinset1 = ax.inset_axes([0.55,0.55,0.4,0.4], facecolor="none")
axinset2 = ax.inset_axes([0.55,0.08,0.4,0.4], facecolor="none")

V_Sr, delta_G, theta_list =  hydrox(T_range)
ax.plot(T_range, V_Sr/x, label="SrO using adsorbed H$_2$O")

V_Sr, delta_G = humid(T_range)
ax.plot(T_range, V_Sr/x, label="SrO using H$_2$O$_{(g)}$ w (2H)'$_{La}$")
axinset1.plot(T_range, V_Sr/x, label="SrO using H$_2$O$_{(g)}$", color = default_colors[1])

V_Sr, delta_G, p_H2 = humid_case1(T_range)
ax.plot(T_range, V_Sr/x, label="SrO using H$_2$O$_{(g)}$ w H$_{2(g)}$")
axinset2.plot(T_range, V_Sr/x, label="SrO using H$_2$O$_{(g)}$ w H$_{2(g)}$", color=default_colors[2])

V_Sr, delta_G, delta_oxygen_list,Ni_H2O_list  = sroh2_case(T_range)
ax.plot(T_range, V_Sr/x, label="Sr(OH)$_{2(g)}$")
axinset2.plot(T_range, V_Sr/x, label="Sr(OH)$_2$", color = default_colors[3])

ax.set_xlim(T_lower_bound, T_upper_bound+1)
ax.set_ylim(0,)
ax.set_xlabel("T[K]")
ax.set_ylabel("$[V\'\'\'_{La}]_{eq}=x_{eq}$ normalised to 1")
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())


axinset1.set_xlim(T_lower_bound, T_upper_bound+1)
axinset1.set_ylim(0,)
axinset1.xaxis.set_minor_locator(AutoMinorLocator())
axinset1.yaxis.set_minor_locator(AutoMinorLocator())
axinset1.yaxis.set_major_formatter(fsf('%.0e'))


axinset2.set_xlim(T_lower_bound, T_upper_bound+1)
axinset2.set_ylim(0,)
axinset2.xaxis.set_minor_locator(AutoMinorLocator())
axinset2.yaxis.set_minor_locator(AutoMinorLocator())
axinset2.yaxis.set_major_formatter(fsf('%.0e'))


ax.legend(loc="upper left", facecolor="none")


fig.savefig(local_path + "figs/humid_cond_summar.png", format="png", dpi=300, transparent=True)

plt.show()
