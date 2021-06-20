import numpy as np

def fDarcy(D, v, rho, mu, eps, f_0=0.1, tol=0.0001):
    """This function allows you calculate the Darcy's friction factor with the Colebrooke-White equation through an iterative method.
    Input your data in consisten units, wether english system or international system, etc.
    
    Esta función te permite calcular el factor de fricción de Darcy con la ecuación de Colebrooke-White a través de un método iterativo.
    Pon tus datos en unidades consistentes, ya sean del sistema internacional o del sistema inglés, etc."""

    Re = (rho * v * D)/mu               # Número de Reynolds
    f = 1/(-2 * np.log10(eps/(3.7*D) + 2.51/(Re*np.sqrt(f_0))))**2       # Valor inicial para primera iteración

    while abs(f - f_0) > tol:
        f_0 = f
        f = 1/(-2 * np.log10(eps/(3.7*D) + 2.51/(Re*np.sqrt(f_0))))**2
    else:
        return f
