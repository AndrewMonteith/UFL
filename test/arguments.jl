using UFL, Test

c = Constant(1.0)

@test ufl_shape(c) === ()
@test c.value === 1.0