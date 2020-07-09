import numpy as np 
from inparams import *
from mesh_refinement import *

"""
Material Routine returns,

1. Material Stiffness Matrix

"""
def material_routine(strain,deltastrain,dt,Q,T,prestress):
    C_e = [[(E/((1+nu)*(1-2*nu)))*(1 - nu), (E/((1+nu)*(1-2*nu)))*nu],[(E/((1+nu)*(1-2*nu)))*nu,(E/((1+nu)*(1-2*nu)))*(1 - nu)]]
    Ct = C_e + ((1 /(1 + (dt/T))) * Q * np.array([2/3, - 1/3, -1/3 , 2/3]).reshape(2,2))
    overstress =  ((1/(1+(dt/T)))*(np.dot(Q*np.array([2/3,-1/3,-1/3,2/3]).reshape(2,2),deltastrain) + prestress))
    stress = np.dot(C_e,strain) + overstress
    return stress,overstress,Ct

