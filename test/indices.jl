using UFL, Test 

i = Index() 
j = Index() 

@test j != i 
@test i < j 
@test i <= j 
@test i <= i

ii = UFL.MultiIndexNode((i,))
@test i ∈ ii 
@test j ∉ ii

ii′ = UFL.MultiIndexNode(())
@test i ∉ ii′

ii′′ = UFL.MultiIndexNode((i, j))
@test all(j ∈ ii′′ for j ∈ (i, j))

@test ii != ii′ 

ij, ji = UFL.MultiIndexNode((i, j)), UFL.MultiIndexNode((j, i))
@test ij != ji

v = [j, i]
sort!(v)
@test v == [i, j]