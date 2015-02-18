def pause():
    garbage = raw_input("[PAUSE]")

from grid import *
from gui import ising_gui
import numpy as np
import matplotlib.pyplot as plt

gui = True

# ---------------------------------------
# Change lattice dimensions here
Ni = 20
Nj = 20
# ---------------------------------------

grid = ising_grid(Ni, Nj)
if gui: display = ising_gui(Ni, Nj)    

M, M_err = [], []
Ch, ChB = [], []
T = np.linspace(1.0, 4.0, 30)

grid.B = 0.00
grid.J = 1.00

# ---------------------------------------
# Change Metropolis iteration scale here
grid.xM = 50
# ---------------------------------------

for T_i in T:
    grid.set_T(T_i)
    M_i, M_std = grid.metropolis()
    if gui: display.draw_grid(grid)
    print "T = %.3e, M = %.3e"%(T_i, M_i)
    M.append(M_i)
    M_err.append(M_std)
    
M = np.array(M)
fig = plt.figure()
ax_1 = fig.add_subplot(111)
ax_1.errorbar(T, M, yerr=M_err, fmt="b.")
plt.xlabel("Temperature (K)")
plt.ylabel("Magnetization")

plt.show()
