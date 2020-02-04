using UFL, Test

i, j = Identity(3), Identity(3)

@test UFL.compute_expr_hash(i) === UFL.compute_expr_hash(j)

s1, s2 = i+j, j+i

# @test UFL.compute_expr_hash(s1) === UFL.compute_expr_hash(s1)
# @test UFL.compute_expr_hash(s1) === UFL.compute_expr_hash(s2)
