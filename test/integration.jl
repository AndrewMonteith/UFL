using UFL, Test


dx = Measure("dx")

@test dx.integral_type === "cell"
@test dx.subdomain_id === "everywhere"

ds = Measure("ds")
@test ds.integral_type === "exterior_facet"
@test ds.subdomain_id === "everywhere"

V = VectorFunctionSpace(UnitSquareMesh(), "CG", 1)

u = UflFunction(V)
c = Constant((0.0, -0.25))

x = det(grad(u))*dx
@test x isa Form 
@test length(x.integrals) === 1

x′= dot(c,u)*dx
@test x′ isa Form 
@test length(x′.integrals) === 1 

@test dx.subdomain_id === "everywhere"
dx′ = dx(4)
@test dx′.subdomain_id === 4

x′′ = dot(c, u)*dx′
@test x′′ isa Form 
@test length(x′′.integrals) === 1


s = x+x′
@test s isa Form
@test length(s.integrals) === 2
@test s.integrals[1] === x.integrals[1]
@test s.integrals[2] === x′.integrals[1]
@test s + 0 === s
@test 0 + s === s