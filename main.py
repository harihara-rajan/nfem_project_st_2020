import numpy as np 
import matplotlib.pyplot as plt 
from elemental_routine import *
from analytical import *
from inparams import *

Q =int(input("Enter the Value of Q: "))
time = np.arange(0,30.1,0.1)
U = np.zeros((num_ele+1,1), dtype=float)
disp_B = []
du = np.zeros((num_ele+1,1), dtype=float)
fext_global = np.zeros((num_ele+1,1), dtype=float)
A = assignment()
over_stress = np.zeros((2,num_ele))
strain = np.zeros((2,num_ele))
load = []
for i in time:
    count = 0
    dt = time[1]- time[0]
    if i <= tl:
        fext_global[0]= (pmax * a * (i))/tl # Load scaling
        load.append((pmax * a * (i))/tl)
    else:
       fext_global[0] = pmax * a
       load.append(pmax * a)
    """
    Newton Raphson Scheme 
    """
    while True: 
        count+=1
        fint_global = np.zeros((num_ele+1,1), dtype=float)
        K_global =np.zeros((num_ele+1,num_ele+1), dtype=float)
        R = np.zeros((num_ele+1,1),dtype=float)
        """
        Calling Element Routine and Assembling Global Matrices
        """
        for i in range(num_ele):
           u = np.array([U[i],U[i+1]]).reshape(2,1)
           du_e = np.array([du[i],du[i+1]]).reshape(2,1)
           prestress = over_stress[:,i].reshape(2,1)
           f_int_element, k_elem, OStress = fint_elements_K_elements(i,u,dt,Q,du_e,prestress)
           over_stress[0,i] = OStress[0]
           over_stress[1,i] = OStress[1]
           fint_global += np.dot(np.transpose(A[i]),f_int_element)
           K_global += A[i].transpose().dot(k_elem.dot(A[i]))
        R = fext_global - fint_global  
        du = np.dot(np.linalg.inv(K_global), R)
        U = U + du
        if np.linalg.norm(R,np.inf) <= 0.005*np.linalg.norm(fint_global,np.inf) or np.linalg.norm(du,np.inf) <= 0.005 * np.linalg.norm(U,np.inf):
            break        
    disp_B.append(U[-1])    
analytical_displacement = analytical()
#print(len(disp_B))
#print(len(time))
plt.plot(time, disp_B)
plt.plot(rnodes, analytical_displacement, c = 'k' ,label = 'Analytical Solution')
plt.scatter(rnodes,U, c='m', marker = 'o', label = 'Numerical Solution')
plt.legend()
plt.show()