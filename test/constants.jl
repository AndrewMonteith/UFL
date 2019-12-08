using UFL, Test


z = Zero(())
@test z.ufl_shape === () 
@test z.ufl_free_indices === () 
@test z.ufl_index_dimensions === () 


z = Zero((3,))
@test z.ufl_shape === (3,)
@test z.ufl_free_indices === () 
@test z.ufl_index_dimensions === ()

i, j = Index(2), Index(4)

z = Zero((), (j, i), Dict(i => 3, j => 5))
@test z.ufl_shape === () 
@test z.ufl_free_indices == (2, 4)
@test z.ufl_index_dimensions === (3, 5)