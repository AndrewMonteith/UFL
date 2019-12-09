A sample implemendation of UFL in Julia

### Current Goals:
Get this code to work:
```julia 
mesh is 2D

    V = VectorFunctionSpace(mesh, "P", 1)

    v = TestFunction(V) <- Shape (2, )

    u = Function(V) <- Shape (2, )

    T = Constant((0, -0.5))

    B = Constant((0, -0.25))

    d = u.geometric_dimension() <- d

    I = Identity(d) (2x2 identity)

    F = I + grad(u) <- (2x2)             # Deformation gradient

    C = F.T*F <- (sum_i F.T[i, j] * F[j, k])                   # Right Cauchy-Green tensor

    # Invariants of deformation tensors

    Ic = tr(C)

    J = det(F)

    # Lamé parameters, "quite squishy"

    mu = Constant(6.3)

    lmbda = Constant(10.0)

    # Stored strain energy density (compressible neo-Hookean model)

    psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2

    # Total potential energy
    Pi = psi*dx # - dot(T, u)*ds(4) - dot(B, u)*dx

    F = derivative(Pi, u, v)
```

### Current Questions for my supervisor:
1. Parameters in the evaluate function. mapping, component & index_values? The evaluate function may not be needed for my side of the project
