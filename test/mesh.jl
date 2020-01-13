using UFL, Test

m = Mesh(triangle)

@test topological_dimension(m) == 2 
@test geometric_dimension(m) == 2