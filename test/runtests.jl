using UFL, Test

@testset "Constant Values" begin include("constantvalue.jl") end
@testset "Algebra" begin include("algebra.jl") end
@testset "Example Form" begin include("exampleform.jl") end