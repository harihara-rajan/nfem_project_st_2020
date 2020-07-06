import numpy as np 
import matplotlib.pyplot as plt 
from elemental_routine import *
from analytical import *
from inparams import *

time = np.arange(0,30.1,0.1)
U = np.zeros((num_ele+1,1))
fext_global = np.zeros((num_ele+1,1))
A = assignment()
for i in time:
    count = 0
    if i <= tl:
        p = ( pmax * a * (i))/tl # Load scaling 
    else:
        p = pmax * a
    fext_global[0] = p
    while True:
        print(i)
        fint_global = np.zeros((num_ele+1,1))
        K_global =np.zeros((num_ele+1,num_ele+1))
        R = np.zeros((num_ele+1,1))
        for i in range(num_ele):
           u = np.array([U[i],U[i+1]]).reshape(2,1) 
           f_int_element, k_elem = fint_elements_K_elements(i,u)
           fint_global += np.dot(np.transpose(A[i]),f_int_element)
           K_global += A[i].transpose().dot(k_elem.dot(A[i])) 
        R = fext_global - fint_global   
        du = np.linalg.solve(K_global,R)
        if np.linalg.norm(R,np.inf) <= 0.005*np.linalg.norm(fint_global,np.inf) or np.linalg.norm(du,np.inf) < 0.005 * np.linalg.norm(U,np.inf):
            break
        U = U + du
        count += 1

# analytical_displacement = analytical()
# plt.plot(rnodes, analytical_displacement, c = 'b' ,label = 'Analytical Solution ')
# plt.scatter(rnodes,U, c='orange', label = 'Numerical Solution')
# plt.legend()
# plt.show()


