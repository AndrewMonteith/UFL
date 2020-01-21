module UFL

include("util.jl")

include("coretypes.jl")
include("indices.jl")
include("uflmacros.jl")

include("iteration.jl")

include("cell.jl")
include("fem.jl")

include("constantvalue.jl")
include("constants.jl")
include("algebra.jl")

include("indexed.jl")
include("differentiation.jl")
# include("indexsum.jl")

include("mesh.jl")
include("functions.jl")
include("arguments.jl")


end # module