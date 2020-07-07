import numpy as np 
from inparams import *
from mesh_refinement import *

"""
Material Routine returns,

1. Material Stiffness Matrix

"""
def elastic_stiffness():
    C_e = [[(E/((1+nu)*(1-2*nu)))*(1 - nu), (E/((1+nu)*(1-2*nu)))*nu],[(E/((1+nu)*(1-2*nu)))*nu,(E/((1+nu)*(1-2*nu)))*(1 - nu)]]
    return np.array(C_e).reshape(2,2)

def material_stiffness(Q,dt):
    C = [[(E/((1+nu)*(1-2*nu)))*(1 - nu), (E/((1+nu)*(1-2*nu)))*nu],[(E/((1+nu)*(1-2*nu)))*nu,(E/((1+nu)*(1-2*nu)))*(1 - nu)]] + ((1 /(1 + dt/T)) * Q * np.array([2/3, - 1/3, -1/3 , 2/3]).reshape(2,2))
    return C

print(material_stiffness(100000,0.5))

