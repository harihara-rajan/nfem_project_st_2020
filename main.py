import numpy as np 
import matplotlib.pyplot as plt 
from elemental_routine import *
from analytical import *
from inparams import *
from mesh_refinement import *
"""
Main program 
------------------------------
takes inputs & results from 
inparams, elemental routine.py, 
,analytical.py & mesh_refinement
------------------------------
Plotting
1. widening of the pipe at the outer radius 'b'
2. convergence study
3. Radial Stress Distribution  
4. Tangential Stress Distribution 
------------------------------
"""
print('For linear elastic case press 0')
print('For visco elastic case press 1')
inputs = int(input("solve for elastic (0) or visco-elatic(1): "))

if inputs == 0:
    Q = 0 
else:
    Q = 100000

time = np.arange(0,30.1,0.1) # time sequence
U = np.zeros((num_ele+1,1), dtype=float) # global nodal displacement
disp_b = [] # array to store displacement at outer radius b, from time = 0 seconds to 30 seconds
disp_a = [] # array to store displacement at inner radius a, from time = 0 seconds to 30 seconds 
du = np.zeros((num_ele+1,1), dtype=float) 
fext_global = np.zeros((num_ele+1,1), dtype=float) # initializing global external force
A = assignment() # calling assignment matrix from the element routine
over_stress = np.zeros((2,num_ele))
stress =  np.zeros((2,num_ele)) # initializing stress array, will be updated 
analytical_displacement = analytical() 

for i in time:
    count = 0
    dt = time[1]- time[0]
    if i <= tl:
        fext_global[0]= (pmax * a * (i))/tl # Load scaling
    else:
       fext_global[0] = pmax * a
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
           u = np.array([U[i],U[i+1]]).reshape(2,1) # initial guess for displacement
           du_e = np.array([du[i],du[i+1]]).reshape(2,1) # 
           prestress = over_stress[:,i].reshape(2,1) # initial guess for over stress at m'th time step
           f_int_element, k_elem, OStress, stress_elemental = fint_elements_K_elements(i,u,dt,Q,du_e,prestress)
           over_stress[0,i] = OStress[0] #copying the overstress of each element to empty over stress array: 
           over_stress[1,i] = OStress[1] 
           stress[0,i] = stress_elemental[0] # 0th row of stress coresponds to radial stress 
           stress[1,i] = stress_elemental[1] # 1st row of stress coresponds to tangential stress
           fint_global += np.dot(np.transpose(A[i]),f_int_element) # global internal force 
           K_global += A[i].transpose().dot(k_elem.dot(A[i])) # global element stiffness matrix
        R = fext_global - fint_global  # finding Residual 
        du = np.dot(np.linalg.inv(K_global), R) # solvin for du
        U = U + du # updating the global displacement 
        """ checkinf for convergence  """
        if np.linalg.norm(R,np.inf) <= 0.005*np.linalg.norm(fint_global,np.inf) or np.linalg.norm(du,np.inf) <= 0.005 * np.linalg.norm(U,np.inf):
            break        
    disp_b.append(U[-1]) # displacements at b: plotting purpose 
    disp_a.append(U[0])  # displacements at a: plotting purpose

"""
Extracting Radial and Tangential stress from the Stress vector(voigt Notation)
"""
radial_stress = stress[0,:]
tangential_stress = stress[1,:]

fig, ax = plt.subplots(1,2)
"""
Convergence Study between Analytical and Numerical solution   - plot  
"""
ax[0].plot(rnodes, analytical_displacement, c = 'k' ,label = 'Analytical Solution')
ax[0].scatter(rnodes,U, c='m', marker = 'o', label = 'Numerical Solution')
ax[0].set(xlabel = 'Radius [mm]', ylabel = 'Displacement [mm]', title = 'Convergence Study')
ax[0].grid()
ax[0].legend()


""" 
Widening of pipe at inner and outer radius - plot 
"""
ax[1].plot(time, disp_a, label = 'Inner radius a')
ax[1].plot(time, disp_b, label = 'Outer radius b')
ax[1].set(xlabel = 'Time [s]', ylabel = 'Displacement [mm]' , title = 'Widening of Pipe at Inner and Outer radius')
ax[1].legend()
ax[1].grid()


fig, ax = plt.subplots(1,2)
"""
Radial Stress Distribution - plot
"""
ax[0].plot(radius, radial_stress, c = 'blue' ,label = 'Radial Stress')
ax[0].set(xlabel = 'Radius [mm]', ylabel = 'Radial Stress [MPa]', title = 'Radial Stress Distribution')
ax[0].legend()
ax[0].grid()
"""
Tangential Stress Distribution - plot 
"""
ax[1].plot(radius, tangential_stress, c = 'orange' ,label = 'Tangential Stress')
ax[1].set(xlabel = 'Radius [mm]', ylabel = 'Tangential Stress [Mpa]', title = 'Tangential Stress Distribution')
ax[1].grid()
ax[1].legend()
plt.show()