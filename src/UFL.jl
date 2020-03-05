module UFL

include("util.jl")

include("coretypes.jl")
include("uflmacros.jl")
include("indices.jl")
include("traversal.jl")

include("cell.jl")
include("fem.jl")

include("geometry.jl")
include("differentiation.jl")

include("constantvalue.jl")
include("tensors.jl")
include("indexing.jl")

include("algebra.jl")
include("tensoralgebra.jl")
include("mathfunctions.jl")

include("mesh.jl")
include("functions.jl")
include("arguments.jl")

include("integration.jl")
include("corefunctions.jl")
include("precedence.jl")
include("apply_differentiation.jl")
include("mapdag.jl")

include("benchmark_differentiation.jl")
end # module