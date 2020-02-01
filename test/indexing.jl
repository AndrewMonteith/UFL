using UFL, Test

v = as_tensor([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

i, j = Index(), Index() 

x = v[i, j]

@test typeof(x) === Indexed
@test ufl_shape(x) === () 
@test ufl_free_indices(x) == (i, j)
@test ufl_index_dimensions(x) === (3, 3)

k = FixedIndex(3)

x2 = v[k, i]

@test typeof(x2) === Indexed 
@test ufl_shape(x2) === () 
@test ufl_free_indices(x2) == (i,)
@test ufl_index_dimensions(x2) == (3,)

x4 = v[i, i]

@test typeof(x4) === IndexSum
@test ufl_shape(x4) === () 
@test ufl_free_indices(x4) === () 
@test ufl_index_dimensions(x4) === () 

x5 = v[i, :]
@test typeof(x5) === ComponentTensor
@test ufl_shape(x5) === (3,)
@test ufl_free_indices(x5) == (i,)
@test ufl_index_dimensions(x5) === (3,)

x6 = v[:, :]
@test x6 === v

x7 = v[i, k]
@test ufl_shape(x7) === () 
@test ufl_free_indices(x7) == (i,)
@test ufl_index_dimensions(x7) === (3,)


data = [[
    [1, 2, 3], [4, 5, 6], [7, 8, 9],
    [1, 2, 3], [4, 5, 6], [7, 8, 9],
    [1, 2, 3], [4, 5, 6], [7, 8, 9],
],
[
    [1, 2, 3], [4, 5, 6], [7, 8, 9],
    [1, 2, 3], [4, 5, 6], [7, 8, 9],
    [1, 2, 3], [4, 5, 6], [7, 8, 9],
],
[
    [1, 2, 3], [4, 5, 6], [7, 8, 9],
    [1, 2, 3], [4, 5, 6], [7, 8, 9],
    [1, 2, 3], [4, 5, 6], [7, 8, 9],
]]

v2 = as_tensor(data)

x8 = v2[i, FixedIndex(3), 3]
@test ufl_shape(x8) === ()
@test ufl_free_indices(x8) == (i,)
@test ufl_index_dimensions(x8) === (3,)