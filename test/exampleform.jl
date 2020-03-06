using Test, UFL 

mesh = UnitTriangleMesh()

V = VectorFunctionSpace(mesh, "CG", 1)

v = TestFunction(V) # <- Shape (2, )

u = UflFunction(V) # <- Shape (2, )

T = Constant((0.0, -0.5))

B = Constant((0.0, -0.25))

d = geometric_dimension(u) # <- d

I = Identity(d) # (2x2 identity)

F = I + grad(u) # <- (2x2)              Deformation gradient

C = F.T*F # <- (sum_i F.T[i, j] * F[j, k])                   # Right Cauchy-Green tensor

# # Invariants of deformation tensors

Ic = tr(C)

J = det(F)

# # Lamé parameters, "quite squishy"

mu = Constant(6.3)

lmbda = Constant(10.0)

# # Stored strain energy density (compressible neo-Hookean model)
psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))^2

# # Total potential energy

Pi = Ic*dx

# Pi = psi*dx #- dot(T, u)*ds(4) - dot(B, u)*dx
# Pi = dot(T, u)*ds(4)

Pi = dot(B, u)*dx - dot(T, u)*ds(4)

F = derivative(Pi, u; du=v)

F_lowered = apply_algebra_lowering(F)

println("Lowered:", F_lowered)
F′ = apply_derivatives(F_lowered)

println("Final:", F′)