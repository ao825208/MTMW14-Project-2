"""
MTMW14: Project 2 - A simplified model of ocean gyres and the Gulf Stream
Student No. 30825208
Parameters
"""

# Parameters
L = 1e6                 # Computational domain of grid
f_0 = 1e-4              # Initial coriolis parameter coefficient
B = 1e-11               # Beta-plane coriolis parameter coefficient
g = 10                  # Gravitation acceleration
gamma = 1e-6            # Linear drag coefficient
rho = 1e3               # Uniform density
H = 1e3                 # Resting depth of fluid (assumed constant)
tau_0 = 0.2             # Wind stress acting on fluid surface

# Dimensions
d = 5e4                 # Grid spacing (both x and y directions)
dt = 175                # Timestep
nd = int(L/d)           # Grid size (No. of grid-points in x, y)
t = 60*60*24            # Time (1 day)
nt = int(t/dt)          # Number of time steps