using Test, UFL 

cell = triangle 

f = FiniteElement("CG"; cell=cell, degree=1)

@test fem_family(f) === "Lagrange"
@test fem_degree(f) === 1
@test (isempty ∘ fem_value_shape)(f)
@test (isempty ∘ fem_ref_value_shape)(f)
@test fem_mapping(f) === "identity"

v = VectorElement("CG"; cell=cell, degree=1)

@test fem_family(v) === "Lagrange" 
@test fem_degree(v) === 1
@test fem_value_shape(v) === (2,)
@test fem_ref_value_shape(v) === (2,)
@test fem_mapping(v) === "identity"

m = MixedElement(f, v)

@test fem_family(m) === "Mixed"
@test fem_degree(m) === 1
@test fem_value_shape(m) === (3,)
@test fem_ref_value_shape(m) === (3,)
@test fem_mapping(m) === "identity"

scalar_f = FiniteElement("Real"; degree = 0)
@test fem_family(scalar_f) === "Real"
@test fem_degree(scalar_f) === 0 
@test (isempty ∘ fem_value_shape)(scalar_f)
@test (isempty ∘ fem_ref_value_shape)(scalar_f)
@test fem_mapping(scalar_f) === "identity"


scalar_v = VectorElement("Real", degree = 0, dim=3)