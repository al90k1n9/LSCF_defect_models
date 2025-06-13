from lib.auxilliary_functions import *
import matplotlib.pyplot as plt
import numpy as np

params = delta_oxygen_interpolater(plot=1)


print("T=973 and the delta is")
print(params[0]*973 + params[1])

