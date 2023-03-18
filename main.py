"""
MTMW14: Project 2 - A simplified model of ocean gyres and the Gulf Stream
Student No. 30825208
Main Code
"""

from numpy import pi, sin, cos, exp, linspace, meshgrid, sqrt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
from parameters import *
from analytical import *
from numerical import *
from energy import *
from numba_exp import *
import time

# Total runtime of code: (approx. 5 mins)

# Task C (Analytic solution)

def analytic():
    """Function that plots the analytic solution derived by Mushgrave (1985)
    onto individual contour plots for each steady state variable."""

    # Grid set-up in both dimensions
    x = linspace(0, L, int(nd))
    y = linspace(0, L, int(nd))
    
    # Defining the solution
    X, Y = meshgrid(x, y)
    u_st, v_st, eta_st = analytic_soln(X, Y)

    # Building a figure
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (15, 4))
    fig.suptitle("Task C: Analytical solution", fontsize = 25)

    # Plotting a contour graph for n_st 
    cp1 = ax1.contourf(X, Y, eta_st)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size = '5%', pad = 0.15)
    ax1.set_title('$\eta$', fontsize = 20)
    ax1.set_xlabel('$x$ [0, L] ($m$)', fontsize = 10)
    ax1.set_ylabel('$y$ [0, L] ($m$)', fontsize = 10)
    fig.colorbar(cp1, cax = cax, orientation = 'vertical',\
                 label = 'Surface elevation ($m$)')

    # Plotting a contour graph for u_st 
    cp2 = ax2.contourf(X, Y, u_st)
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size = '5%', pad = 0.15)
    ax2.set_title('$u$', fontsize = 20)
    ax2.set_xlabel('$x$ [0, L] ($m$)', fontsize = 10)
    ax2.set_ylabel('$y$ [0, L] ($m$)', fontsize = 10)
    fig.colorbar(cp2, cax = cax, orientation = 'vertical', \
                 label = 'Horizontal velocity ($ms^{-1}$)')

    # Plotting a contour graph for v_st 
    cp3 = ax3.contourf(X, Y, v_st)
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size = '5%', pad = 0.15)
    ax3.set_title('$v$', fontsize = 20)
    ax3.set_xlabel('$x$ [0, L] ($m$)', fontsize = 10)
    ax3.set_ylabel('$y$ [0, L] ($m$)', fontsize = 10)
    fig.colorbar(cp3, cax = cax, orientation = 'vertical', \
                 label = 'Vertical velocity ($ms^{-1}$)')

    fig.tight_layout()
    
analytic()
print("Task C completed.")

# Task D (Forward-backward time scheme)

def numerical(t):
    """Function that defines the plots for the forward-backward time scheme to
    solve the shallow water momentum equations derived by Matsuno (1966)."""
    
    # Values of u, v, eta
    uNew, vNew, etaNew, total_energy_time = FBT_scheme(d, t, dt)
    
    # Building a figure
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize = (25, 20))
    fig.suptitle("Task D: Forward-backward time scheme model", fontsize = 35)
    
    # Plotting 1D arrays for the x and y axis
    x = np.linspace(0, L, int(nd + 2))
    y = np.linspace(0, L, int(nd + 2))
    
    # Plot 1: u vs. x
    ax1.plot(x, uNew[3, :])
    ax1.set_xlabel('$x$ [0, L] ($m$)', fontsize = 20)
    ax1.set_ylabel('$u$ ($ms^{-1}$)', fontsize = 20)
    ax1.tick_params(axis = 'both', which = 'major', labelsize = 20)
    ax1.tick_params(axis = 'both', which = 'minor', labelsize = 20)
    ax1.set_title('$u$ vs $x$', fontsize = 30)
    
    # Plot 2: v vs. y
    ax2.plot(y, vNew[:, 4])
    ax2.set_xlabel('$y$ [0, L] ($m$)', fontsize = 20)
    ax2.set_ylabel('$v$ ($ms^{-1}$)', fontsize = 20)
    ax2.tick_params(axis = 'both', which = 'major', labelsize = 20)
    ax2.tick_params(axis = 'both', which = 'minor', labelsize = 20)
    ax2.set_title('$v$ vs $y$', fontsize = 30)
    
    # Plot 3: eta vs. x
    gyre_mid = int((nd + 2)/2)
    ax3.plot(x, etaNew[gyre_mid, :])
    ax3.set_xlabel('$\eta$ ($m$)', fontsize = 20)
    ax3.set_ylabel('$y$ [0, L] ($m$)', fontsize = 20)
    ax3.tick_params(axis = 'both', which = 'major', labelsize = 20)
    ax3.tick_params(axis = 'both', which = 'minor', labelsize = 20)
    ax3.set_title('$\eta$ vs $x$', fontsize = 30)
    
    # Plot 4: 2D Contour plot of eta
    cp = ax4.contourf(x, y, etaNew[:, :])
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes('right', size = '5%', pad = 0.15)
    ax4.set_title('$\eta$', fontsize = 30)
    ax4.set_xlabel('$x$ [0, L] ($m$)', fontsize = 20)
    ax4.set_ylabel('$y$ [0, L] ($m$)', fontsize = 20)
    ax4.tick_params(axis = 'both', which = 'major', labelsize = 20)
    ax4.tick_params(axis = 'both', which = 'minor', labelsize = 20)
    cbar = fig.colorbar(cp, cax = cax, orientation = 'vertical')
    cbar.ax.tick_params(labelsize = 15)
    cbar.set_label('Elevation ($m$)', fontsize = 20)
    
    fig.tight_layout()
    
    # Prints the value required for eta_0 in the analytical condition
    print("eta_0 is: ", etaNew[int((nd+2)/2), 1])

# D.1 (Sanity check)
numerical(t)       
print("Task D.1 completed.")

# D.2 (Steady state)

# Calculating processing time for numerical scheme (required for task G)
start1 = time.time()
numerical(t*40)
end1 = time.time()
print("Elapsed time for Task D computation = %s" % (end1 - start1),"s")
print("Task D.2 completed.")

# D.3 (Difference plots)

def diff_plots():
    """Function that plots the differences between the numerical and
    the steady state analytical solutions."""
    
    # Defining x and y spaces
    x1 = linspace(0, L, int(nd + 2))
    y1 = linspace(0, L, int(nd + 2))
    
    X1, Y1 = meshgrid(x1, y1)
    
    # Defining solutions from the steady state and numerical models
    u_st, v_st, eta_st = analytic_soln(X1, Y1)
    uNew, vNew, etaNew, total_energy_time = FBT_scheme(d, t*40, dt)
    
    # Difference between the numerical and analytic solutions
    u_diff = uNew - u_st
    v_diff = vNew - v_st
    eta_diff = etaNew - eta_st
    
    # Defining the different energy differences
    energy_diff = total_energy(u_diff, v_diff, eta_diff)
    energy_st = total_energy(u_st, v_st, eta_st)
    energy_num = total_energy(uNew, vNew, etaNew)
    
    # Building a figure
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (15, 4))
    fig.suptitle("Task D.3: Energy differences", fontsize = 25)

    # Plotting a contour graph for eta_diff 
    cp1 = ax1.contourf(x1, y1, eta_diff)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size = '5%', pad = 0.15)
    ax1.set_title('$\eta$', fontsize = 20)
    ax1.set_xlabel('$x$ [0, L] ($m$)', fontsize = 10)
    ax1.set_ylabel('$y$ [0, L] ($m$)', fontsize = 10)
    fig.colorbar(cp1, cax = cax, orientation = 'vertical',\
                 label = 'Surface elevation ($m$)')

    # Plotting a contour graph for u_diff 
    cp2 = ax2.contourf(x1, y1, u_diff)
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size = '5%', pad = 0.15)
    ax2.set_title('$u$', fontsize = 20)
    ax2.set_xlabel('$x$ [0, L] ($m$)', fontsize = 10)
    ax2.set_ylabel('$y$ [0, L] ($m$)', fontsize = 10)
    fig.colorbar(cp2, cax = cax, orientation = 'vertical', \
                 label = 'Horizontal velocity ($ms^{-1}$)')

    # Plotting a contour graph for v_diff 
    cp3 = ax3.contourf(x1, y1, v_diff)
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size = '5%', pad = 0.15)
    ax3.set_title('$v$', fontsize = 20)
    ax3.set_xlabel('$x$ [0, L] ($m$)', fontsize = 10)
    ax3.set_ylabel('$y$ [0, L] ($m$)', fontsize = 10)
    fig.colorbar(cp3, cax = cax, orientation = 'vertical', \
                 label = 'Vertical velocity ($ms^{-1}$)')
    
    fig.tight_layout()

diff_plots()
print("Task D.3 completed.")

# Task E (Energy method)

def energy_method(d, dt):
    """Function that plots a time series graph of the evolution of energy from
    the numerical model."""
    
    # Defining x and y spaces
    x1 = linspace(0, L, int(nd + 2))
    y1 = linspace(0, L, int(nd + 2))
    
    X1, Y1 = meshgrid(x1, y1)
    
    # Defining solutions from the steady state and numerical models
    u_st, v_st, eta_st = analytic_soln(X1, Y1)
    uNew, vNew, etaNew, total_energy_time = FBT_scheme(d, t*40, dt)
    
    # Difference between the numerical and analytic solutions
    u_diff = uNew - u_st
    v_diff = vNew - v_st
    eta_diff = etaNew - eta_st
    
    # Defining the different energy differences
    energy_diff = total_energy(u_diff, v_diff, eta_diff)
    energy_st = total_energy(u_st, v_st, eta_st)
    energy_num = total_energy(uNew, vNew, etaNew)

    # Plotting the graph of total energy vs. time
    time_axis = linspace(0, t*40, int(t*40/dt))
    
    # Building a figure
    fig, ax = plt.subplots(figsize = (8, 5))
    fig.suptitle("Task E: Total Energy time series", fontsize = 25)
    
    # Plotting the solution
    ax.plot(time_axis, total_energy_time)
    ax.set_ylabel('Energy ($J$)')
    ax.set_xlabel('Time ($s$)')
    ax.set_title('Total energy vs. time')
    
    fig.tight_layout()
    
    # Determining the error measurement between the steady state and numerical
    # solutions
    print("Energy error measurement =", energy_diff, "J")
    # Determining the total energy amount when the model reaches steady state
    print("Steady state total energy value = ", total_energy_time[-1], "J")
    
energy_method(d, dt)

# Halving the grid spacing and reducing the time step to satisfy 2D CFL
energy_method(d*0.5, 100)
print("Task E completed.")

# Task G.2 (Optimisation of computational aspects with numba)

# Filters numba warnings
import warnings
warnings.filterwarnings('ignore')

# Calculating processing time for numerical scheme with numba
start2 = time.time()
FBT_schemeG(d, t*40, dt);
end2 = time.time()
print("Elapsed time for Task G computation = %s" % (end2 - start2), "s")
print("Overall time difference =", (end2 - start2) - (end1 - start1), "s")
print("Task G.2 completed.")