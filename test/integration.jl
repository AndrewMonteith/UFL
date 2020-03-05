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

s = det(grad(u))
s2 = dot(c, u)

x = s*dx
x′= s2*dx