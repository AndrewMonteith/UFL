using BenchmarkTools

export run_benchmark

function build_random_tree(N)
    i1 = Identity(3) + Identity(3)
    r = rand()
    i2 = as_tensor([[r, r, r], [r, r, r], [r, r, r]]) + Identity(3)

    rootnode = Identity(3) + Identity(3)
    
    for i in 1:N
        r = rand()
        if r < -0.5 
            rootnode = rootnode + i1 
        elseif r < 0
            rootnode = i2 + rootnode 
        elseif r < 0.5 
            rootnode = rootnode + Identity(3)
        else
            rootnode = Identity(3) + rootnode 
        end 
    end

    rootnode
end

function count_nodes(x::AbstractExpr, operands::VarTuple{Int})
    if isempty(operands)
        1
    else
        1 + length(operands)
    end
end


do_benchmark_1(tree) = map_expr_dag(tree, count_nodes)
do_benchmark_2(tree) = map_expr_dag(tree, count_nodes, Int)

function run_benchmark()
    suite = BenchmarkGroup()

    n = 400
    println("building with $(n) nodes")

    suite["building_tree"] = @benchmarkable build_random_tree($n)
    # suite["without_type"] = @benchmarkable do_benchmark_1(x) setup=(x=build_random_tree($n))
    # suite["with_type"] = @benchmarkable do_benchmark_2(x) setup=(x=build_random_tree($n))

    tune!(suite)

    BenchmarkTools.run(suite, verbose=true, seconds=10)
end
