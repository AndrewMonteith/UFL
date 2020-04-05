from ufl.corealg.map_dag import map_expr_dag
from ufl.corealg.multifunction import MultiFunction
from ufl.corealg.traversal import pre_traversal, post_traversal
from ufl.algorithms.apply_algebra_lowering import apply_algebra_lowering
from ufl.algorithms.apply_derivatives import apply_derivatives


from firedrake import * 
import random

mesh = UnitSquareMesh(2, 2)
s = SpatialCoordinate(mesh)
(x, y) = s 

V = VectorFunctionSpace(mesh, "CG", 1)

u = TrialFunction(V)
v = TestFunction(V)

# def count_nodes(root):
#     s = 0
#     for x in pre_traversal(root):
#         s += 1
#     return s

# def count_nodes_post(root):
#     s = 0 
#     for x in post_traversal(root):
#         s += 1 
#     return s

def build_random_tree(n, x, u, v):
    atom1 = tr(grad(u))
    atom2 = dot(u, v)
    atom3 = x**2 
    atom4 = ln(x**2) 

    def random_bool():
        return random.random() <= 0.5

    root = atom3

    for i in range(round((n-5)/21)):
        for _ in range(4):
            r = random.random()
            if r < 0.25:
                if random_bool():
                    root = root + atom1 
                else:
                    root = atom1 + root 
            elif r < 0.5:
                if random_bool():
                    root = root * atom2 
                else:
                    root = atom2 * root 
            elif r < 0.75:
                if random_bool():
                    root = root + atom3 
                else:
                    root = atom3 + root 
            else:
                if random_bool():
                    root = root + atom4
                else:
                    root = atom4 + root 
    
    return root

tree_100 = apply_algebra_lowering(build_random_tree(100, x, u, v))
tree_200 = apply_algebra_lowering(build_random_tree(200, x, u, v))
tree_300 = apply_algebra_lowering(build_random_tree(300, x, u, v))
tree_500 = apply_algebra_lowering(build_random_tree(500, x, u, v))
tree_1000 = apply_algebra_lowering(build_random_tree(1000, x, u, v))
tree_5000 = apply_algebra_lowering(build_random_tree(5000, x, u, v))
tree_10000 = apply_algebra_lowering(build_random_tree(10000, x, u, v))
tree_20000 = apply_algebra_lowering(build_random_tree(20000, x, u, v))