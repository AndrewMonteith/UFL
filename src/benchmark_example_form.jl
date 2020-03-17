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
psi = (Ic - 3)*(mu/2) - mu*ln(J) + (lmbda/2)*(ln(J))^2

# Total potential energy
dx = Measure("dx")
ds = Measure("ds")

# F_lowered = apply_algebra_lowering(F)
# println("Lowered:", F_lowered)
# F′ = apply_derivatives(F_lowered)

# println("Final:", F′)
function benchmark(f::Form)
    apply_derivatives(apply_algebra_lowering(f))
end

using BenchmarkTools

function do_benchmarks()
    suite = BenchmarkGroup()
    Pi = psi*dx - dot(T, u)*ds(4) - dot(B, u)*dx
    F = derivative(Pi, u; du=v)

    suite["form"] = @benchmarkable benchmark($F) # Julia:~433ms Python:271ms

    tune!(suite)

    BenchmarkTools.run(suite, verbose=true, seconds=10)
end