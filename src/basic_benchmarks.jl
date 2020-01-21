include("UFL.jl")

using Main.UFL

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

function random_tree_benchmark()
    x = build_random_tree()
    print("done")
end

function do_benchmarks(randomtree)
    # build up random expression tree

    s = 0
    for node in pre_order_traversal(randomtree)
        s += 1
    end

    s
end

function do_benchmarks_(randomtree)
    s = 0 

    pre_order_traversal_((_) -> s += 1, randomtree) 

    s
end


for i in 1:200
    println("Benchark $(i)")
    s = build_random_tree(1000000)

    @time do_benchmarks(s)
    @time do_benchmarks_(s)
end