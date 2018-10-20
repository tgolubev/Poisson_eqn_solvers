# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 21:44:03 2018

@author: Tim

This solves the 1D Poisson equation by finite differences, assuming a constant
dielectric constant. A tridiagonal matrix is set up and solved.Con

Inputs: number of mesh points, left and right voltage boundary conditions,
dielectric constant.
"""

import numpy as np, time
import matplotlib
import matplotlib.pyplot as plt  # allows to make matlab style plots!

def thomasSolve(matrix, rhs):
    
    x = np.zeros(num_elements+2)  # array for the solution
    
   # Forward substitution
    for i in range(2, matrix.num_elements + 1): # recall range uses [ ) ]
        cdiag_ratio = matrix.lower_diag[i-1]/matrix.main_diag[i-1]
        matrix.main_diag[i] -= cdiag_ratio * matrix.upper_diag[i-1]
        rhs[i] -= cdiag_ratio * rhs[i-1]
        
    # Backward substitution
    x[num_elements] = rhs[num_elements]/matrix.main_diag[num_elements] #lin eqn. corresponding to last row
    for i in range(num_elements, 1, -1):  # 3rd argument is the step (iterate down)
        x[i-1] = (rhs[i-1] - x[i]*matrix.upper_diag[i-1])/matrix.main_diag[i-1]
        
    return x
    
        
class Poisson_matrix():
    
    def __init__(self, epsilon, num_elements):
        self.epsilon = epsilon  # save epsilon as attribute of the class
        self.num_elements = num_elements
        self.main_diag = -2*epsilon*np.ones(num_elements+1)
        self.upper_diag = epsilon*np.ones(num_elements)
        self.lower_diag = self.upper_diag  #poisson matrix is symmetric
        
        

# System setup
        
#num_elements =  int(input("Enter number of mesh points ")) 
#V_left_BC = float(input("Voltage at the left boundary "))
#V_right_BC = float(input("Voltage at the right boundary "))
#epsilon = float(input("Relative dielectric constant "))
        
num_elements = 99
V_left_BC = 0
V_right_BC = 0.9
epsilon = 3.8

start = time.time()

poisson_mat = Poisson_matrix(epsilon, num_elements)  # initialize object

# setup right hand side--> for now just hard code in a dipole, later will 
# read the charges from file
rhs =  np.zeros(num_elements+1)  # make it +1 b/c I want to index from 1...
rhs[24] = 0.681
rhs[25] = -0.681

# apply boundary conditions
rhs[1] -= epsilon*V_left_BC
rhs[num_elements] -= epsilon*V_right_BC

V = thomasSolve(poisson_mat, rhs)

# add the boundary values (which didn't need to be solved for)
V[0] = V_left_BC
V[num_elements+1] = V_right_BC

end = time.time()
print(V)
print(f"CPU time: {end - start}")
# C++ version of this with 9999 elements runs in ~0.0032 sec

# make plot
fig, ax = plt.subplots()  # make a figure
ax.plot(V)
ax.set(xlabel = 'Position', ylabel='Voltage (V)', title = 'Poisson solver result')
ax.grid()  # add a grid
fig.savefig("poisson_result.jpg")
plt.show()
