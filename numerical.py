"""
MTMW14: Project 2 - A simplified model of ocean gyres and the Gulf Stream
Student No. 30825208
Numerical Model (Forward-backward time scheme with Arakawa C-grid)
"""

from numpy import *
from parameters import *
from energy import *

def FBT_scheme(d, t, dt):
    """Forward-backward time scheme that solves for new values of eta, u, v the 
    shallow water equations, also using the Arakawa C grid to define the
    derivatives dudx, dvdy and detady."""
    
    # Number of time steps
    nt = int(t/dt)          

    # Creating arrays for new values of u, v and eta
    uNew = np.zeros((nd + 2, nd + 2))
    vNew = np.zeros((nd + 2, nd + 2))
    etaNew = np.zeros((nd + 2, nd + 2))
    
    # Creating arrays for old values of u, v and eta
    uOld = np.zeros((nd + 2, nd + 2))
    vOld = np.zeros((nd + 2, nd + 2))
    etaOld = np.zeros((nd + 2, nd + 2))
    
    total_energy_time = np.zeros([nt])
    
    # For each odd time step, u is calculated before v
    Odd_nt = True
    
    # Stepping through time
    for n in range(nt): 
        # Stepping through y direction
        for j in range(1, nd + 1): 
            # Stepping through x direction
            for i in range(1, nd + 1): 
                # Arakawa C-grid implementation
                dudx = (uOld[j, i + 1] - uOld[j, i])/d
                dvdy = (vOld[j + 1, i] - vOld[j, i])/d 
                etaNew[j, i] = etaOld[j, i] - (H*dt)*(dudx + dvdy) 
        # Odd timestep
        if Odd_nt:  
            # Calculating new u
            # Stepping through y direction
            for j in range(1, nd + 1): 
                # Wind stress calculated at each j step
                tau = [-np.cos(np.pi*j*d/L)*tau_0, 0] 
                # Stepping through x direction
                for i in range(1, nd + 1):              
                    detadx = (etaNew[j, i] - etaNew[j, i - 1])/d
                    fv = (f_0 + B*j*d)*(vOld[j, i] + vOld[j + 1, i] +\
                                          vOld[j, i - 1] +\
                                              vOld[j + 1, i - 1])/4
                    uNew[j, i] = uOld[j, i] +  dt*fv  - g*dt*detadx -\
                        gamma*dt*uOld[j, i] + tau[0]*dt/(rho*H)
            # Calculating new v
            # Stepping through y direction
            for j in range(1, nd + 1): 
                # Wind stress calculated at each j step
                tau = [-np.cos(np.pi*j*d/L)*tau_0 , 0] 
                # Stepping through x direction
                for i in range(1, nd + 1): 
                    detady = (etaNew[j, i] - etaNew[j - 1, i])/d
                    fu = (f_0 + B*j*d)*(uNew[j, i] + uNew[j, i + 1] +\
                                          uNew[j - 1, i] +\
                                              uNew[j - 1, i + 1])/4
                    vNew[j, i] = vOld[j, i] - dt*fu  - g*dt*detady -\
                        gamma*dt*vOld[j, i] + tau[1]*dt/(rho*H)
                        
        # For each even time step, v is calculated before u
        else: 
            # Calculating new v
            # Stepping through y direction
            for j in range(1, nd + 1):
                # Wind stress calculated at each j step
                tau = [-np.cos(np.pi*j*d/L)*tau_0, 0] 
                # Stepping through x direction  
                for i in range(1, nd + 1):              
                    detady = (etaNew[j, i] - etaNew[j - 1, i])/d
                    fu = (f_0 + B*j*d)*(uOld[j, i] + uOld[j, i + 1] +\
                                          uOld[j - 1, i] +\
                                              uOld[j - 1, i + 1])/4
                    vNew[j, i] = vOld[j, i] - dt*fu  - g*dt*detady -\
                        gamma*dt*vOld[j, i] + tau[1]*dt/(rho*H)
            # Calculating new u
            # Stepping through y direction
            for j in range(1, nd + 1):
                # Wind stress calculated at each j step
                tau = [-np.cos(np.pi*j*d/L)*tau_0, 0] 
                # Stepping through x direction 
                for i in range(1, nd + 1):              
                    detadx = (etaNew[j, i] - etaNew[j, i - 1])/d
                    fv = (f_0 + B*j*d)*(vNew[j, i] + vNew[j + 1, i] +\
                                          vNew[j, i - 1] +\
                                              vNew[j + 1, i - 1])/4
                    uNew[j, i] = uOld[j, i] +  dt*fv  - g*dt*detadx -\
                        gamma*dt*uOld[j, i] + tau[0]*dt/(rho*H)
        
        # Calculate energy integral
        total_energy_time[n] = total_energy(uNew, vNew, etaNew)
        
        # Toggling between odd/even
        Odd_nt = not Odd_nt 
        
        # Resetting variables
        uOld = uNew
        vOld = vNew
        etaOld = etaNew
    
    return uNew, vNew, etaNew, total_energy_time