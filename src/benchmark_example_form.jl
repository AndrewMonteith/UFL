using BenchmarkTools

function build_form()
    mesh = UnitTriangleMesh()

    V = VectorFunctionSpace(mesh, "CG", 1)

    v = TestFunction(V) # <- Shape (2, )

    u = UflFunction(V) # <- Shape (2, )

    T = Constant((0.0, -0.5))

    B = Constant((0.0, -0.25))

    d = geometric_dimension(u) # <- d

    I = Identity(d) # (2x2 identity)

    F = I + grad(u) # <- (2x2)              Deformation gradient

    C = trans(F)*F # <- (sum_i F.T[i, j] * F[j, k])                   # Right Cauchy-Green tensor

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
    
    Pi = psi*dx - dot(T, u)*ds(4) - dot(B, u)*dx
    F = derivative(Pi, u; du=v)

    F
end


function do_benchmarks()
    suite = BenchmarkGroup()
    # # Stored strain energy density (compressible neo-Hookean model)

    # Total potential energy

    F = build_form()
    F_lowered = apply_algebra_lowering(F)
    F′ = apply_derivatives(F_lowered)

    println("--- Doing Pullback")
    F_pulled = apply_function_pullback(F′)
    F′′ = apply_derivatives(F_pulled)
    

    # suite["lowering"] = @benchmarkable apply_algebra_lowering($F)
    # suite["differentiation_1"] = @benchmarkable apply_derivatives($F_lowered)
    suite["pulled"] = @benchmarkable apply_function_pullback($F′)
    # suite["differentiation_2"] = @benchmarkable apply_derivatives($F_pulled)

    tune!(suite)

    BenchmarkTools.run(suite, verbose=true, seconds=5)
end