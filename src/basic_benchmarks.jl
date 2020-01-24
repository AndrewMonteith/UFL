using BenchmarkTools

function build_random_tree(N)
    rootnode = Identity(3) + Identity(3)
    
    for i in 1:N
        if randn() < 0 
            rootnode = rootnode + Identity(3)
        else
            rootnode = Identity(3) + rootnode 
        end
    end

    rootnode
end

function random_tree_benchmark(N)
    x = build_random_tree(N)
    nothing
end

function count_nodes(expr::AbstractExpr)
    to_visit::Vector{AbstractExpr} = [expr] 
    s = 0

    while !isempty(to_visit)
        e = pop!(to_visit)
        append!(to_visit, ufl_operands(e))

        s += 1
    end

    s
end

function do_benchmarks_1(tree::AbstractExpr)
    s = 0

    pre_order_traversal_(tree) do x
        s += 1
    end

    s
end

function do_benchmarks_2(tree::AbstractExpr)
    s = 0

    pre_order_traversal_(tree) do x
        s += 1
    end

    s
end

function do_benchmarks_3(tree::AbstractExpr) 
    s = 0 

    @pre_order_traversal_m(tree, begin
        s += 1 
    end)

    s 
end

function do_benchmarks_4(tree::AbstractExpr) 
    s = 0 

    @pre_order_traversal_m2(tree, begin 
        s += 1 
    end)

    s 
end

function do_benchmarks_5(tree)
    s = 0
    for node in Main.UFL.pre_order_traversal(tree) 
        s += 1
    end
    s
end

function run_benchmark()
    suite = BenchmarkGroup()

    println("building with 10_000 nodes")

    suite["building-array"] = @benchmarkable build_random_tree(10_000)
    suite["function-inbuilt-array"] = @benchmarkable do_benchmarks_1(x) setup=(x=build_random_tree(10_000))
    suite["function-capacity-array"] = @benchmarkable do_benchmarks_2(x) setup=(x=build_random_tree(10_000))
    suite["macro-inbuilt-array"] = @benchmarkable do_benchmarks_3(x) setup=(x=build_random_tree(10_000))
    suite["macro-capacity-array"] = @benchmarkable do_benchmarks_4(x) setup=(x=build_random_tree(10_000))
    suite["iterator-inbuilt-array"] = @benchmarkable do_benchmarks_5(x) setup=(x=build_random_tree(10_000))
    suite["raw"] = @benchmarkable count_nodes(x) setup=(x=build_random_tree(10_000))

    tune!(suite)

    BenchmarkTools.run(suite, verbose=true, seconds=5)
end