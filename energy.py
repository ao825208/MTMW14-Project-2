"""
MTMW14: Project 2 - A simplified model of ocean gyres and the Gulf Stream
Student No. 30825208
The Energy method
"""

import numpy as np
from parameters import *

def total_energy(u, v, eta):
    """Function that numerically calculates the total energy of the
    perturbation from the resting ocean."""
    
    # Integration arrays
    x = np.linspace(0, L, len(u[1]))
    y = np.linspace(0, L, len(u[0]))
    
    # Energy method equation
    def E(u, v, eta):
        return 0.5*rho*(H*(u**2 + v**2) + g*eta**2)
    
    E = E(u, v, eta)
    
    # Integrating energy equation in both x and y arrays
    total_energy = np.sum(E)*d*d
    
    return total_energy