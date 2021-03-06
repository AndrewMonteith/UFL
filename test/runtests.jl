using UFL, Test

@testset "Indices" begin include("indices.jl") end
@testset "Constant Values" begin include("constantvalue.jl") end
@testset "Algebra" begin include("algebra.jl") end
@testset "Meshs" begin include("mesh.jl") end
@testset "FEM" begin include("fem.jl") end
@testset "Misc" begin include("misc.jl") end
@testset "Hashing" begin include("hashing.jl") end
@testset "Iteration" begin include("traversal.jl") end
@testset "Tensors" begin include("tensors.jl") end
@testset "Indexing" begin include("indexing.jl") end
@testset "Map Dag" begin include("mapdag.jl") end
@testset "Differentiation" begin include("differentiation.jl") end
@testset "Tensor Algebra" begin include("tensoralgebra.jl") end
@testset "Integration" begin include("integration.jl") end
@testset "Algebraic Lowering" begin include("algebra_lowering.jl") end
@testset "Function Pullback" begin include("function_pullback.jl") end
@testset "Example Form" begin include("exampleform.jl") end