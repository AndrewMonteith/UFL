using BenchmarkTools

function build_random_tree(N::Int)
    root = Identity(3) + Identity(3)

    for i ∈ 1:N 
        if rand() < -0.5 
            root = root + Identity(3)
        elseif rand() < 0
            root = Zero((3, 3)) * root 
        elseif rand() < 0.5 
            root = root + as_tensor([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        elseif rand() < 1 
            root = (-Identity(3)) * root
        end
    end 

    root
end 

function benchmark_1(tree)
    s = 0
    i = 0
    @UFL.pre_order_traversal for x ∈ tree
        i += 1
        s += length(ufl_shape(x))
    end 
    (i, s)
end


function run_benchmarks()
    suite = BenchmarkGroup()

    n = 500
    println("building with $(n) nodes")

    suite["building-tree"] = @benchmarkable build_random_tree($n)
    suite["trait-iterators"] = @benchmarkable benchmark_1(x) setup = (x = build_random_tree($n))

    tune!(suite)

    BenchmarkTools.run(suite, verbose = true, seconds = 10)
end