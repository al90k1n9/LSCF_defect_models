import numpy as np

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
