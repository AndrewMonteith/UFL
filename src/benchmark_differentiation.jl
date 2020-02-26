using BenchmarkTools

function build_random_tree(N::Int)
    s = SpatialCoordinate(3, 3)

    (x, y, z) = s 

    root = x + x^2 

    for i ∈ 1:N 
        term = x^(rand(1:100000))
        if rand() < 0 
            root = root + term
        else
            root = term + root
        end 
    end 

    root
end 

function benchmark_1(tree::AbstractExpr)
    apply_derivatives(tree)
end

function run_benchmarks()
    suite = BenchmarkGroup()

    n = 10000
    println("building with $(n) nodes")

    suite["building-tree"] = @benchmarkable build_random_tree($n) # Julia:~433ms Python:271ms
    suite["differentiation"] = @benchmarkable benchmark_1(x) setup=(x=build_random_tree($n)) # Julia: ~6.2ms Python: ~220ms

    tune!(suite)

    BenchmarkTools.run(suite, verbose=true, seconds=10)
end