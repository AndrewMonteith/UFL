using UFL, Test

@testset "Constant Values" begin include("constantvalue.jl") end
@testset "Algebra" begin include("algebra.jl") end
@testset "Meshs" begin include("mesh.jl") end
@testset "FEM" begin include("fem.jl") end 
@testset "Constants" begin include("arguments.jl") end
@testset "Misc" begin include("misc.jl") end 
# @testset "Iteration Tests" begin include("iteration.jl") end
@testset "Example Form" begin include("exampleform.jl") end