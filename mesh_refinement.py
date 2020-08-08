a = 50 #inner radius 
b = 100 # outer radius 
num_ele = 10 #int(input("Number of elements: ")) # num of elements 

##### Mesh Refinement Factor ######
mrf = 2
if num_ele == 1:
    q = mrf
else:
    q = mrf**(1/(num_ele - 1))

dr=(b-a)*((1-q)/(1-mrf*q))
rnode = a
rnodes = [a]

for i in range(num_ele):
    rnode = rnode + dr 
    rnodes.append(rnode)
    dr = dr * q

"""
This is for plotting tangential and radial stress distribution with the radius
radius used in main program for plotting purpose 
"""
radius = []
for i in range(num_ele):
    radius.append((rnodes[i] + rnodes[i+1])/ 2)

