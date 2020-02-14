using UFL, Test

i, j = Identity(3), Identity(3)

@test i.ufl_hash_code === j.ufl_hash_code
@test hash(i) === hash(j)

s1, s2 = i+j, j+i

@test s1.ufl_hash_code === s2.ufl_hash_code
@test hash(s1) === hash(s2)

exprs = Set{UFL.AbstractExpr}()
push!(exprs, i)
push!(exprs, s1)

@test i ∈ exprs 
@test j ∈ exprs 

@test s1 ∈ exprs 
@test s2 ∈ exprs 

k = Identity(4)
@test k ∉ exprs
