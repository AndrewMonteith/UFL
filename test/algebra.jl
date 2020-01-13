using Test, UFL 


i = Identity(3)
j = Identity(3)

@assert ufl_shape(i) == (3, 3)
@assert ufl_shape(j) == (3, 3)

i_sum = i + j 

@test ufl_shape(i_sum) == (3, 3)