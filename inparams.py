"""
------------------------------
Given input parameters 
------------------------------
"""

E = 200000 # --> youngs Modulus
nu = 0.20  # --> poisons ratio
T =  3     # --> time like parameter for load scaling 
a = 50     # --> inner radius of the pipe
b = 100    # --> outer radius of the pipe 
pmax = 140 # --> max pressure [MPa]
tl = 6     # --> time in second in which load has to reach 140 MPa
ts = 30    # --> after 6 seconds the load shoud be constant
