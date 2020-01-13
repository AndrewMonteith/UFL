using Test, UFL 

cell = triangle 

f = FiniteElement("CG", cell, 1)

@test fem_family(f) === "Lagrange"
@test fem_degree(f) === 1
@test fem_value_shape(f) === ()

v = VectorElement("CG", cell, 1)

@test fem_family(v) === "Lagrange" 
@test fem_degree(v) === 1
@test fem_value_shape(v) === (2,)

m = MixedElement(f, v)

@test fem_family(m) === "Mixed"
@test fem_degree(m) === 1
@test fem_value_shape(m) === (3,)