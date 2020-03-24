using UFL, Test

V = VectorFunctionSpace(UnitSquareMesh(), "CG", 1)
u = UflFunction(V)
c = Constant((0, -0.25))

x = det(grad(u))# *dx
lowered_x = apply_algebra_lowering(x)
r = r"\(grad\(w\_\d\)\)\[0, 0\] \* \(grad\(w\_\d\)\)\[1, 1\] \+ -1 \* \(grad\(w\_\d\)\)\[0, 1\] \* \(grad\(w_\d+\)\)\[1, 0\]"
@test occursin(r, string(lowered_x))
@test (isempty ∘ ufl_free_indices)(lowered_x)

x = dot(c, u)# *dx 
lowered_x = apply_algebra_lowering(x)
r = r"sum\_\{i\_\d+\} w\_\d+\[i\_\d+\] \* w\_\d+\[i\_\d+\]"
@test occursin(r, string(lowered_x))
@test (isempty ∘ ufl_free_indices)(lowered_x)

x = trans(grad(u))
lowered_x = apply_algebra_lowering(x)
r = r"\{ A \| A\_\{i\_\d+, i\_\d+\} = \(grad\(w\_\d+\)\)\[i\_\d+, i\_\d+\] \}"
@test (isempty ∘ ufl_free_indices)(lowered_x) && (isempty ∘ ufl_index_dimensions)(lowered_x)
@test ufl_shape(lowered_x) === (2, 2)
@test occursin(r, string(lowered_x))
@test (isempty ∘ ufl_free_indices)(lowered_x)

x = tr(grad(u))
lowered_x = apply_algebra_lowering(x)
r = r"sum\_\{i\_\d+\} \(grad\(w\_\d+\)\)\[i\_\d+, i\_\d+\]" 
@test occursin(r, string(lowered_x))
@test (isempty ∘ ufl_free_indices)(lowered_x)

x = grad(tr(grad(u)))
lowered_x = apply_algebra_lowering(x)
r = r"grad\(sum\_\{i\_\d+\} \(grad\(w\_\d+\)\)\[i\_\d+, i\_\d+\]\)"
@test occursin(r, string(lowered_x))
@test (isempty ∘ ufl_free_indices)(lowered_x)