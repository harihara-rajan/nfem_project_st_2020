import numpy as np 
from inparams import *
from mesh_refinement import *

"""
Material Routine returns,

1. Tangent Stiffness Matrix

"""

Q = 0
dt = 0.01
T = 1
def C():
    G = E / (1+nu)*(1-2*nu)
    C = G * np.array([1-nu,nu,nu,1-nu]).reshape(2,2)
    return C


C = [[(E/((1+nu)*(1-2*nu)))*(1 - nu), (E/((1+nu)*(1-2*nu)))*nu],[(E/((1+nu)*(1-2*nu)))*nu,(E/((1+nu)*(1-2*nu)))*(1 - nu)]] + (1 /(1 + dt/T)) * Q * np.array([2/3, - 1/3, -1/3 , 2/3]).reshape(2,2)