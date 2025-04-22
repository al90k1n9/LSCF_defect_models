import numpy as np

acell = 14.641598078 #Bohr
bohr2angstrom = 0.52177

red2cart = np.asarray([1,1,5]) * acell * bohr2angstrom

og_o1 = np.asarray([7.4076004134E-01, 2.1591012494E-01, 4.1481175529E-01])#58
og_h1 = np.asarray([8.1539519077E-01, 3.1418832058E-01, 4.1970138000E-01])#87

#indices  1 are bonded and indices 2 are bonded

og_o2 = np.asarray([7.5741977865E-01, 1.9588505177E-01, 9.8492142701E-01])#82
og_h2 = np.asarray([8.2533191524E-01, 2.9543639395E-01, 9.7782933190E-01])#88

#gives you the direction in which your hydrogen will be, with oxygen atoms fixed
vect1 = np.multiply(-og_o1+og_h1, red2cart)
vect2 = np.multiply(-og_o2+og_h2, red2cart)

og_dist1 = np.sqrt(np.sum(np.multiply(vect1, red2cart)**2))
og_dist2 = np.sqrt(np.sum(np.multiply(vect2, red2cart)**2))

unit_vect1 = vect1/og_dist1
unit_vect2 = vect2/og_dist2

u_list = [-0.15, 0.15, 0.1, -0.1, 0.05, -0.05, 0.025, -0.025]
for u in u_list:
    new_h1 = np.multiply(og_o1, red2cart) + (1+u) * og_dist1 * unit_vect1
    new_h2 = np.multiply(og_o2, red2cart) + (1+u) * og_dist2 * unit_vect2
    
    
    
    print(np.multiply(new_h1, 1/red2cart), " new reduced position 87 ", u)
    print(np.multiply(new_h2, 1/red2cart), " new reduced position 88 ", u)
    print("\n \n")

