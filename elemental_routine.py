import numpy as np 
from mesh_refinement import *
from material_routine import *

def N_mat():
    chi = 0
    N=[(1/2) * (1 - chi), (1/2) * (1 + chi)]
    return N 


def jacob(rnodes, r1, r2):
    return (r2 - r1) / (2)


def B_mat(rnodes, r1, r2):
    B = [[-1 / (r2 - r1), +1 / (r2 - r1)],[+1 / (r1 + r2),+1 / (r1 + r2)]]
    return B


# def ele_stiff():
#     """Function for calculating element stiffness matrix"""
#     B = B_mat(rnodes, rnodes[i], rnodes[i+1])
#     J = jacob(rnodes, rnodes[i], rnodes[i+1])
#     N = N_mat()
#     r = [rnodes[i], rnodes[i+1]]
#     k1 = np.transpose(B).dot(C.dot (B))
#     k2 = J * np.dot(np.transpose(N),np.transpose(r))
#     ke = 2*k1*k2
#     return ke

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

def fint_elements_K_elements(i,U):
    B = B_mat(rnodes, rnodes[i], rnodes[i+1])
    J = jacob(rnodes, rnodes[i], rnodes[i+1])
    N = N_mat()
    r = [rnodes[i], rnodes[i+1]]
    strain = (np.dot(B,U))
    stress = np.dot(C,strain)
    s1 = np.dot(np.transpose(B), stress)
    s2 = J * np.dot(np.transpose(N),np.transpose(r))
    fint_element = 2*s1*s2
    ###########
    k1 = np.transpose(B).dot(C.dot (B))
    k2 = J * np.dot(np.transpose(N),np.transpose(r)) 
    k_element =2*k1*k2 
    return fint_element,k_element 

    



