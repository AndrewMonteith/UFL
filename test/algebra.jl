using Test, UFL 


i = Identity(3)
j = Identity(3)

i_sum = i + j 

@test ufl_shape(i_sum) == (3, 3)