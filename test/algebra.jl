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

@test (isempty ∘ ufl_operands)(s)
@test (isempty ∘ ufl_shape)(s)

x1 = as_tensor([[1, 2, 3]])
x2 = as_tensor([[1], [2], [3]])

m = x1*x2
@test ufl_shape(m) === (1, 1)
@test (isempty ∘ ufl_free_indices)(m)
@test (isempty ∘ ufl_index_dimensions)(m)

m = 3*x1
@test ufl_shape(m) === (1, 3)
@test (isempty ∘ ufl_free_indices)(m)
@test (isempty ∘ ufl_index_dimensions)(m)

m = -x1 
@test m isa ComponentTensor 
@test (isempty ∘ ufl_free_indices)(m)
@test (isempty ∘ ufl_index_dimensions)(m)