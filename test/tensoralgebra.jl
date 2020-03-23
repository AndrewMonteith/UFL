using UFL, Test 

i = Identity(3)

j = trans(i)

@test ufl_shape(j) === (3, 3)
@test (isempty ∘ ufl_free_indices)(j)