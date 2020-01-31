using UFL, Test

v = as_tensor([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

i, j = Index(), Index() 

x = v[i, j]

@test typeof(x) === Indexed
@test ufl_shape(x) === () 
@test ufl_free_indices(x) == (i, j)
@test ufl_index_dimensions(x) == (3, 3)

k = FixedIndex(3)

x2 = v[k, i]

@test typeof(x2) === Indexed 
@test ufl_shape(x2) === () 
@test ufl_free_indices(x2) == (i,)
@test ufl_index_dimensions(x2) == (3,)