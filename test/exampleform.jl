using Test, UFL 

# VectorFunctionSpace does't exist
V = FunctionSpace()

v = TestFunction(V)

# nor does Function 
u = TrialFunction(V)

@test ufl_shape(v) === (2,)
@test ufl_shape(u) === (2,)

# The Constant class in UFL takes (domain, shape)
# The Constant in the example snippet provides neither?

# T = Constant((0, -0.5))
# B = Constant((0, -0.25))

d = geometric_dimension(u)

@test d === 2

I = Identity(d)

