using UFL, Test

i, j = Index(), Index()
@test i.id === 0
@test j.id === 1

z = Zero(())
@test z.ufl_shape === () 
@test z.ufl_free_indices === () 
@test z.ufl_index_dimensions === () 
@test z == 0

z = Zero((3,))
@test z.ufl_shape === (3,)
@test z.ufl_free_indices === () 
@test z.ufl_index_dimensions === ()
@test z == Zero((3, ))


i, j = Index(2), Index(4)

z = Zero((), (j, i), Dict(i => 3, j => 5))
@test z.ufl_shape === () 
@test z.ufl_free_indices == (2, 4)
@test z.ufl_index_dimensions === (3, 5)


i = Identity(3)
j, k = FixedIndex(3), FixedIndex(2)
@test i.ufl_shape === (3, 3)
@test i[3, 3] == IntValue(1)
@test i[3, 2] == IntValue(0)
@test i[j, j] == IntValue(1)
@test i[j, k] == IntValue(0)