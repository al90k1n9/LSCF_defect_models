from lib.SrOH2_models import * 
import matplotlib.pyplot as plt

x = 0.4 #molar fraction of Sr
x_O2 = 0.21
x_H2O = 0.08
P = 1 #atm
T_lower_bound = 700 #K
T_upper_bound = 1400 #K
T_range = np.arange(T_lower_bound, T_upper_bound)

p_O2 = x_O2 * P
p_H2O = x_H2O * P

fig, ax = plt.subplots(layout="constrained")

V_Sr = case1(T_range)
ax.plot(T_range, V_Sr, label="case1")

V_Sr = case2(T_range)
ax.plot(T_range, V_Sr, label="case2")

V_Sr = case3(T_range)
ax.plot(T_range, V_Sr, label="case3")

V_Sr = case4(T_range)
ax.plot(T_range, V_Sr, label="case4")


V_Sr = case5(T_range)
ax.plot(T_range, V_Sr, label="case5")

#print(case2(1000))
#print(case3(1000))
#print(case4(1000))
#print(case5(1000))
ax.set_title("Sr(OH)$_2$ formation")
ax.set_xlabel("T[K]")
ax.set_ylabel("[V\'\'\'$_{La}$]")


ax.set_xlim(left=T_lower_bound,right=T_upper_bound)
ax.set_ylim(0,)
ax.legend(loc="upper left")


ax.legend()
plt.show()