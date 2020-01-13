module UFL

include("coretypes.jl")
include("indices.jl")
include("ufltype.jl")


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


end # module