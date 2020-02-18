from ufl import * 
from ufl.corealg.map_dag import map_expr_dag
from ufl.corealg.multifunction import MultiFunction

import random

def build_random_tree(N):     
    i1 = Identity(3) + Identity(3)
    r = random.random()
    i2 = as_tensor([[r, r, r], [r, r, r], [r, r, r]]) + Identity(3)

    rootnode = Identity(3) + Identity(3)
    
    for i in range(N):
        r = random.random()
        if r < 0.25:
            rootnode = rootnode + i1 
        elif r < 0.5:
            rootnode = i2 + rootnode 
        elif r < 0.75:
            rootnode = rootnode + Identity(3)
        else:
            rootnode = Identity(3) + rootnode 

    return rootnode
