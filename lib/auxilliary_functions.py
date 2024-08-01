import numpy as np
R = 8.314 #J.K.mol-1

def pH2_giver(T,p_H2O, p_O2):
    #delta mu calculation
    T_0 = 298 #K
    p_0 = 1 #bar


    cp_O2 = 29.4 #K J/mol
    delta_hf298_O2 = 0 #kJ.mol-1
    s0_O2 = 205.2 #J.mol-1.K


    cp_H2 = 28.8 #K J/mol
    delta_hf298_H2 = 0 #kJ.mol-1
    s0_H2O = 130.680 #J.mol-1.K

    cp_H2O = 33.6 #K J/mol
    delta_hf298_H2O = -241.826 #kJ.mol-1
    s0_H2O = 188.835 #J.mol-1.K-1

    t=T/1000
    h_mat = np.array([t,  t**2/2, t**3/3, t**4/4, -1/t, 1, 0, -1])
    s_mat = np.array([np.log(t), t,  t**2/2.0,  t**3/3, -1./(2 * t**2), 0, 1, 0])
    if T<1000 and T>=700:
        shomate_coeffs = [[30.09200, 6.832514, 6.793435, -2.534480,0.082139, -250.8810, 223.3967, -241.8264],
                            [33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797, 172.707974, 0.0],
                            [30.03235, 8.772972, -3.988133, 0.788313, -0.741599, -11.32468, 236.1663, 0.0]]
        shomate_coeffs = np.asarray(shomate_coeffs)

    elif T>=1000 and T<1700:
        shomate_coeffs = [[30.09200, 6.832514, 6.793435, -2.534480,0.082139, -250.8810, 223.3967, -241.8264],
                            [18.563083, 12.257357, -2.859786, 0.268238, 1.977990, -1.147438, 156.288133, 0.0],
                            [30.03235, 8.772972, -3.988133, 0.788313, -0.741599, -11.32468, 236.1663, 0.0]]
        shomate_coeffs = np.asarray(shomate_coeffs)
    else:
        return "temperature out of range"
    h = np.dot(shomate_coeffs, h_mat * 1000) #J/mol
    s = np.dot(shomate_coeffs, s_mat) #J K /mol

    stoichiometry = np.asarray([1,-1,-1/2])
    h = h * stoichiometry
    s = s * stoichiometry
    

    hf = np.sum(h) + delta_hf298_H2O * 1000
    gf = hf - T * np.sum(s)

    ph2=p_H2O/np.sqrt(p_O2)*np.exp(gf/(R*T))
    return ph2

def cube_root_complex(z):
    #python has problems taking powers of complex functions. so the necessary function, here
    norm = (z.real**2 + z.imag**2)**(1/6)
    phase = np.arctan(z.imag/z.real)/3
    #print(phase, "phase")
    return complex(norm*np.cos(phase), norm*np.sin(phase))

def cubic_model(a,b,c,d):
    #FUNCTION TO FIND CUBIC ROOTS
    determinant = 18 *a*b*c*d - 4*b**3*d + b**2 * c**2 - 4*a*c**3 - 27*a**2 * d**2
    #print(N, "N")
    #print(determinant, "determinant") #If the determinant < 0, then 1 real + 2 complex root, if determinant > 0, then 3 distinc real roots, if = 0 then multiple roots

    delta_zero = b**2 - 3*a*c
    delta_one = 2*b**3 - 9*a*b*c + 27*a**2*d 
    inner_sqrt = delta_one**2 - 4 * delta_zero**3
    re_inner_sqrt = inner_sqrt
    if inner_sqrt >= 0:
        c_minus = np.cbrt((delta_one - np.sqrt(delta_one**2 - 4 * delta_zero**3))/2)
        c_plus = np.cbrt((delta_one + np.sqrt(delta_one**2 - 4 * delta_zero**3))/2)
    else:
        inner_sqrt = np.abs(inner_sqrt)
        c_minus = cube_root_complex(delta_one/2 - complex(0,np.sqrt(inner_sqrt)/2))
        c_plus = cube_root_complex(delta_one/2 - complex(0,np.sqrt(inner_sqrt)/2))

    xi = complex(-1/2, np.sqrt(3)/2)
    xi_coeffs = np.asarray([1,xi, 1/xi])
    xi_powers = np.asarray([xi_coeffs**0, xi_coeffs, xi_coeffs**2])
    terms_minus = -1/(3*a)*np.asarray([b, c_minus, delta_zero/c_minus])
    terms_plus = -1/(3*a)*np.asarray([b, c_plus, delta_zero/c_plus])
    solutions_plus =np.dot( xi_powers, terms_plus)
    solutions_minus = np.dot(xi_powers, terms_minus)
    #print(solutions_minus[0] -solutions_plus[0], "difference in solutions ", solutions_minus[0])
    return solutions_plus
