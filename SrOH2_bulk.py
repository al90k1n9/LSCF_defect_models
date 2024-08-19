from lib.SrOH2_models import * 


x = 0.4 #molar fraction of Sr
x_O2 = 0.21
x_H2O = 0.08
P = 1 #atm

p_O2 = x_O2 * P
p_H2O = x_H2O * P


print(case1(1000))
print(case2(1000))
print(case3(1000))
print(case4(1000))