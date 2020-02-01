using Test, UFL 


i = Identity(3)
j = Identity(3)

@assert ufl_shape(i) == (3, 3)
@assert ufl_shape(j) == (3, 3)

i_sum = i + j 

@test ufl_operands(i_sum) == (i, j)
@test ufl_shape(i_sum) == (3, 3)
@test i_sum isa Sum 

s = ScalarValue(1)

@test ufl_operands(s) === ()
@test ufl_shape(s) === () 