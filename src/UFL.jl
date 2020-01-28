module UFL

include("util.jl")

include("coretypes.jl")
include("indices.jl")
include("uflmacros.jl")

include("iteration.jl")


include("cell.jl")
include("fem.jl")

include("constantvalue.jl")
include("tensors.jl")
include("indexing.jl")

# include("constants.jl") Never used anywhere
include("algebra.jl")

# include("differentiation.jl") Not implemented yet

include("mesh.jl")
include("functions.jl")
include("arguments.jl")

include("basic_benchmarks.jl")

end # module