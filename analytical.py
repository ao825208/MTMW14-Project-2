"""
MTMW14: Project 2 - A simplified model of ocean gyres and the Gulf Stream
Student No. 30825208
Analytical Solution
"""

from numpy import pi, sin, cos, exp, linspace, meshgrid, sqrt
from parameters import *

def analytic_soln(x, y):
    """Function that determines the analytical solution for the shallow water
    equations derived by Mushgrave (1985)."""
    
    # Constant of integration determined from Task D
    eta_0 = 0.07583681292138715
    
    # Constants required to detemine solutions
    epsilon = gamma/(L*B)
    a = (-1 - (sqrt(1 + (2*pi*epsilon)**2)))/(2*epsilon)
    b = (-1 + (sqrt(1 + (2*pi*epsilon)**2)))/(2*epsilon)
    
    # Functions required to determine solutions
    def f1(x):
        return pi*(1 + ((exp(a) - 1)*(exp(b*x))\
                       + (1 - exp(b))*exp(a*x))/(exp(b) - exp(a)))
            
    def f2(x):
        return ((exp(a) - 1)*b*(exp(b*x)) +\
                (1 - exp(b))*a*exp(a*x))/(exp(b) - exp(a))
            
    # Determining values for the steady state solutions of each parameter along
    # the x and y axis
    u_st = (-tau_0/(pi*gamma*rho*H))*(f1(x/L))*cos(pi*y/L)
    v_st = (tau_0/(pi*gamma*rho*H))*(f2(x/L))*sin(pi*y/L)
    eta_st = eta_0 + (tau_0/(pi*gamma*rho*H))*((f_0*L)/g)*\
        ((gamma/(f_0*pi))*(f2(x/L))*cos(pi*y/L) +\
         (1/pi)*(f1(x/L))*sin(pi*y/L)*(1 + (B*y/f_0)\
                           + (B*L/(f_0*pi))*cos(pi*y/L)))
            
    return u_st, v_st, eta_st