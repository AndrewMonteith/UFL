using Test, UFL 

cell = triangle 

f = FiniteElement("CG"; cell=cell, degree=1)

@test fem_family(f) === "Lagrange"
@test fem_degree(f) === 1
@test (isempty ∘ fem_value_shape)(f)

v = VectorElement("CG"; cell=cell, degree=1)

@test fem_family(v) === "Lagrange" 
@test fem_degree(v) === 1
@test fem_value_shape(v) === (2,)

m = MixedElement(f, v)

@test fem_family(m) === "Mixed"
@test fem_degree(m) === 1
@test fem_value_shape(m) === (3,)

scalar_f = FiniteElement("Real"; degree = 0)
@test fem_family(scalar_f) === "Real"
@test fem_degree(scalar_f) === 0 
@test (isempty ∘ fem_value_shape)(scalar_f)


scalar_v = VectorElement("Real", degree = 0, dim=3)