export do_benchmark

using BenchmarkTools 

function build_random_tree(N::Int)
    root = as_ufl(1)

    for i ∈ 1:N 
        x = rand()
        if x < -0.5
            root = root + tr(Identity(3))
        elseif x < 0 
            root = tr(Identity(3)) + root 
        elseif x < 0.5 
            root = root + tr(Identity(3)+Identity(3))
        else 
            root = tr(Identity(3)+Identity(3)) + root 
        end 
    end

    root
end

function benchmark_1(expr::AbstractExpr)
    apply_algebra_lowering(expr)
end


function do_benchmark() 
    suite = BenchmarkGroup()

    n = 10
    println("building with $n nodes")

    suite["building-tree"] = @benchmarkable build_random_tree($n) # Julia:~433ms Python:271ms
    suite["lowering"] = @benchmarkable benchmark_1(x) setup=(x=build_random_tree($n)) # Julia: ~6.2ms Python: ~220ms

    tune!(suite)

    BenchmarkTools.run(suite, verbose=true, seconds=10)
end