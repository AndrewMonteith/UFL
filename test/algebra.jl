using Test, UFL 


i = Identity(3)
j = Identity(3)

@assert ufl_shape(i) == (3, 3)
@assert ufl_shape(j) == (3, 3)

i_sum = i + j 
@test ufl_shape(i_sum) == (3, 3)
@test i_sum isa Sum 


# Subtraction is defined via -1 * <node>
# i_sub = i - j
# @test ufl_shape(i_sub) == (3, 3)