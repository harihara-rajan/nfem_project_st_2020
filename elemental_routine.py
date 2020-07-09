import numpy as np 
from mesh_refinement import *
from material_routine import *
from inparams import *

def N_mat():
    chi = 0
    N=[(1/2) * (1 - chi), (1/2) * (1 + chi)]
    return N 

def jacob(rnodes, r1, r2):
    return (r2 - r1) / (2)

def B_mat(rnodes, r1, r2):
    B = [[-1 / (r2 - r1), +1 / (r2 - r1)],[+1 / (r1 + r2),+1 / (r1 + r2)]]
    return B

def assignment():
    """Function for calculating Assignment Matrix"""
    num_nod = num_ele + 1
    a = np.zeros((2,num_nod))
    A = []  
    for i in range(0,num_ele):
        a = np.zeros((2,num_nod))
        a[0][i] = 1
        a[1][i+1] = 1
        A.append(a)
    return A

def fext_element(p):
    f_ext_ele = np.zeros((num_ele+1,1))
    f_ext_ele[0,0] = -1 * p * a
    return f_ext_ele

def fint_elements_K_elements(i,u,dt,Q,du_e,prestress): #,du_e,prestress,dt
    B = B_mat(rnodes, rnodes[i], rnodes[i+1])
    J = jacob(rnodes, rnodes[i], rnodes[i+1])
    N = N_mat()
    r = [rnodes[i], rnodes[i+1]]
    strain = (np.dot(B,u)) # elemental strain
    deltastrain = (np.dot(B,du_e))
    stress,over_stress,Ct=material_routine(strain,deltastrain,dt,Q,T,prestress)
    #stress = np.dot(elastic_stiffness(),strain) + over_stress(dt,Q,deltastrain,prestress)
    # print(stress)
    """B.transpose*Stress*N.transpose*r.transpose*weight(2)*J - for calculating elemental internal force(following three lines)"""
    s1 = np.dot(np.transpose(B), stress)
    s2 = J * np.dot(np.transpose(N),np.transpose(r))
    fint_element = 2*s1*s2
    """ B.transpose*C*B*J*N.transpose*r.transpose - for calculating the element stiffness matrix (following three lines)"""
    k1 = np.transpose(B).dot(Ct.dot (B))
    k2 = J * np.dot(np.transpose(N),np.transpose(r)) 
    k_element =2*k1*k2 
    return fint_element, k_element,over_stress




