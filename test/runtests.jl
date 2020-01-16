using UFL, Test

@testset "Constant Values" begin include("constantvalue.jl") end
@testset "Algebra" begin include("algebra.jl") end
@testset "Mesh testing" begin include("mesh.jl") end
@testset "FEM Testing" begin include("fem.jl") end 
@testset "Constant" begin include("arguments.jl") end
@testset "Example Form" begin include("exampleform.jl") end
@testset "Misc Tests" begin include("misc.jl") end 