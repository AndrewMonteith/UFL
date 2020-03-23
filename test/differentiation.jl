using UFL, Test 

s = SpatialCoordinate(3, 3)

(x, y, z) = s 

@test (x, y, z) isa Tuple{Indexed, Indexed, Indexed}

# Whilst we cannot test for the equality of some of the indicies in the strings, if the general form of the derivative 
# strings are similar enough then we can consider them equal

# -- Algebraic Derivatives 

diff = apply_derivatives(grad(x))
expected = r"^\{ A \| A_\{i\_\d+\} = I\[1, i\_\d+\] \}"
@test occursin(expected, string(diff))

diff = apply_derivatives(grad(x + x))
expected = r"^\(\{ A \| A\_\{i\_\d+\} = I\[1, i\_\d+\] }\) \+ \(\{ A \| A\_\{i_\d+\} = I\[1, i\_\d+\] \}\)"
@test occursin(expected, string(diff))

diff = apply_derivatives(grad(2*x))
expected = r"\{ A \| A\_\{i\_\d+\} = 2 \* \(\{ A \| A\_\{i\_\d+\} = I\[1, i\_\d+\] \}\)\[i\_\d+\] \}"
@test occursin(expected, string(diff))

diff = apply_derivatives(grad(x^2))
expected = r"\{ A \| A\_\{i\_\d+\} = x\[\d] \* \(\{ A \| A\_\{i\_\d+\} = 2 \* \(\{ A \| A\_\{i\_\d+\} = I\[1, i\_\d+\] \}\)\[i\_\d+\] \}\)\[i\_\d+\] \}"
@test occursin(expected, string(diff))

# -- Differentiation with Indicies 
V = VectorFunctionSpace(UnitSquareMesh(), "CG", 1)
u = UflFunction(V)
c = Constant((0, -0.25))

i = Index() 
x = grad(u)[i, i]
diff = apply_derivatives(x)
r = r"sum\_\{i\_\d+\} \(grad\(w\_\d+\)\)\[i\_\d+, i\_\d+\]"
@test occursin(r, string(diff))

x = grad(u)
diff = apply_derivatives(x)
@test occursin(r"grad\(w\_\d+\)", string(diff))

x = grad((tr ∘ grad)(u)/2)
diff = apply_derivatives(apply_algebra_lowering(x))
r=r"\{ A \| A\_\{i\_\d+\} = \(sum\_\{i\_\d+\} \(\{ A \| A\_\{i\_\d+\} = \(grad\(grad\(w\_\d+\)\)\)\[i\_\d+, i\_\d+, i\_\d+\] \}\)\)\[i\_\d+\] \/ 2 \}"
@test occursin(r, string(diff))