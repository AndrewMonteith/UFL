using UFL, Test

str_integrand(f, n) = string(f.integrals[n].integrand)

V = VectorFunctionSpace(UnitSquareMesh(), "CG", 1)
u = UflFunction(V)
c = Constant((0, -0.25))

x = det(grad(u))*dx
lowered_x = apply_algebra_lowering(x)
r = r"\(grad\(w\_\d\)\)\[0, 0\] \* \(grad\(w\_\d\)\)\[1, 1\] \+ -1 \* \(grad\(w\_\d\)\)\[0, 1\] \* \(grad\(w_\d+\)\)\[1, 0\]"
@test occursin(r, str_integrand(lowered_x, 1))

x = dot(c, u)*dx 
lowered_x = apply_algebra_lowering(x)
r = r"sum\_\{i\_\d+\} w\_\d\[i\_\d+\] \* w\_\d\[i\_\d+\]"
@test occursin(r, str_integrand(lowered_x, 1))