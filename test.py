from ufl.corealg.map_dag import map_expr_dag
from ufl.corealg.multifunction import MultiFunction

from firedrake import * 
import random

mesh = UnitCubeMesh(5, 5, 5)
s = SpatialCoordinate(mesh)
(x, y, z) = s 

def build_random_tree(N, x):
    root = x + x**2 

    for i in range(N): 
        term = x**(random.randint(1,100000))
        if random.random() < 0.5:
            root = root + term
        else:
            root = term + root

    return root

# def build_random_tree(N):     
#     i1 = Identity(3) + Identity(3)
#     r = random.random()
#     i2 = as_tensor([[r, r, r], [r, r, r], [r, r, r]]) + Identity(3)

#     rootnode = Identity(3) + Identity(3)
    
#     for i in range(N):
#         r = random.random()
#         if r < 0.25:
#             rootnode = rootnode + i1 
#         elif r < 0.5:
#             rootnode = i2 + rootnode 
#         elif r < 0.75:
#             rootnode = rootnode + Identity(3)
#         else:
#             rootnode = Identity(3) + rootnode 

#     return rootnode
