using UFL, Test 

s = SpatialCoordinate(3, 3)

(x, y, z) = s 

@test (x, y, z) isa Tuple{Indexed, Indexed, Indexed}

# Whilst we cannot test for the equality of some of the indicies in the strings, if the general form of the derivative 
# strings are similar enough then we can consider them equal

diff = apply_derivatives(grad(x))
@test occursin(r"^\{ A \| A_\{i\_\d+\} = I\[1, i\_\d+\] \}", string(diff))

diff = apply_derivatives(grad(x + x))
@test occursin(r"^\(\{ A \| A\_\{i\_\d+\} = I\[1, i\_\d+\] }\) \+ \(\{ A \| A\_\{i_\d+\} = I\[1, i\_\d+\] \}\)", string(diff))

diff = apply_derivatives(grad(2*x))
@test occursin(r"\{ A \| A\_\{i\_\d+\} = 2 \* \(\{ A \| A\_\{i_\d+\} = I\[1, i\_22\] \}\)\[i\_\d+\] \}", string(diff)) === true

diff = apply_derivatives(grad(x^2))
@test occursin(r"\{ A \| A\_\{i\_\d+\} = x\[\d] \* \(\{ A \| A\_\{i\_\d+\} = 2 \* \(\{ A \| A\_\{i\_\d+\} = I\[1, i\_\d+\] \}\)\[i\_\d+\] \}\)\[i\_\d+\] \}", 
               string(diff))