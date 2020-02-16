module UFL

include("util.jl")

include("coretypes.jl")
include("uflmacros.jl")
include("indices.jl")
include("traversal.jl")

include("cell.jl")
include("fem.jl")

include("constantvalue.jl")
include("tensors.jl")
include("indexing.jl")

include("algebra.jl")

include("mesh.jl")
include("functions.jl")
include("arguments.jl")

include("corefunctions.jl")
include("mapdag.jl")

include("benchmark_accessor_patterns.jl")
end # module