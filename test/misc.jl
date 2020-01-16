using Test, UFL 

i = Identity(3)
@test as_ufl(i) === i 

s = ScalarValue(3)
@test as_ufl(3) == s


s2 = ScalarValue(convert(Float64, 4))
@test as_ufl(convert(Float64, 4)) == s2