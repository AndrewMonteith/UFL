module UFL

include("util.jl")

include("coretypes.jl")
include("uflmacros.jl")
include("indices.jl")
include("traversal.jl")

include("cell.jl")
include("fem.jl")

include("mesh.jl")
include("geometry.jl")

include("exprcontainers.jl")

include("constantvalue.jl")
include("indexing.jl")
include("tensors.jl")

include("algebra.jl")
include("tensoralgebra.jl")
include("mathfunctions.jl")

include("functions.jl")
include("arguments.jl")
include("integration.jl")
include("differentiation.jl")

include("precedence.jl")
include("apply_differentiation.jl")
include("mapdag.jl")
include("compound_expressions.jl")
include("apply_algebra_lowering.jl")
include("apply_function_pullback.jl")

include("corefunctions.jl")
include("benchmark_example_form.jl")
end # module