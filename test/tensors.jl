using Test, UFL 

v = as_tensor([1,2,3])
@test ufl_shape(v) === (3,)

v2 = as_tensor([[1,2,3], [4,5,6], [7,8,9]])
@test ufl_shape(v2) === (3, 3)