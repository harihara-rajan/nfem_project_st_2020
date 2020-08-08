import numpy as np 
from inparams import *
from mesh_refinement import *
import matplotlib.pyplot as plt 
"""
analytical 
------------------------------
returns 
analytical displacements 
        u(r) = (1 + nu) * ((pmax * a**2)/(E*(b**2-a**2))) * ((1-(2*nu))*rnodes[i] + (b**2/rnodes[i]))

"""
U = []
def analytical():
    for i in range(num_ele+1):
        u = (1 + nu) * ((pmax * a**2)/(E*(b**2-a**2))) * ((1-(2*nu))*rnodes[i] + (b**2/rnodes[i]))
        U.append(u)
    return U