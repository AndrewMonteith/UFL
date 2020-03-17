using UFL, Test

V = VectorFunctionSpace(UnitSquareMesh(), "CG", 1)

u, v = UflFunction(V), TestFunction(V) 

x = dot(grad(u), grad(v))
x_pulled = apply_function_pullback(x)
r = r"\(grad\(reference\_value\(w\_\d+\)\)\)\s*\.\s*\(grad\(reference_value\(v\_\d+\)\)\)"
@test occursin(r, string(x_pulled))