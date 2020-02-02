using Test, UFL

i, j = Index(), Index()

z = Zero(())
@test (isempty ∘ ufl_shape)(z)
@test (isempty ∘ ufl_free_indices)(z)
@test (isempty ∘ ufl_index_dimensions)(z)
@test z == 0

z = Zero((3,))
@test ufl_shape(z) === (3,)
@test (isempty ∘ ufl_free_indices)(z)
@test (isempty ∘ ufl_index_dimensions)(z) 
@test z == Zero((3, ))


# i, j = Index(2), Index(4)

# z = Zero((), (j, i), Dict(i => 3, j => 5))
# @test z.ufl_shape  
# @test z.ufl_free_indices == (2, 4)
# @test z.ufl_index_dimensions === (3, 5)


# i = Identity(3)
# j, k = FixedIndex(3), FixedIndex(2)
# @test i.ufl_shape === (3, 3)
# @test i[3, 3] == ScalarValue(1::Int)
# @test i[3, 2] == ScalarValue(0::Int)
# @test i[j, j] == ScalarValue(1::Int)
# @test i[j, k] == ScalarValue(0::Int)