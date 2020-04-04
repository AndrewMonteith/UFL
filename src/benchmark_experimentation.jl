# struct NodeCounter <: AbstractMapper 
#     base::BaseMapper 
#     NodeCounter() = new(BaseMapper{Int}())
# end

# (nc::NodeCounter)(t::Terminal) = 1
# (nc::NodeCounter)(op::Operator) = 1 + sum(nc[child] for child ∈ ufl_operands(op))

# count_nodes(x::AbstractExpr) = map_expr_dag(NodeCounter(), x)

random_bool() = rand() <= 0.5

function build_static_type_tree(n::Int, x::AbstractExpr, u::AbstractExpr, v::AbstractExpr)
    # 8 Types
    atom1 = tr(grad(u)) # 3 nodes
    atom2 = dot(u, v) # 3 nodes
    atom3 = x^2 # 5 nodes
    atom4 = ln(x^2) # 6 nodes

    root = atom3

    for i ∈ 1:(round((n-5)/21))
        for _ ∈ 1:4
            r = rand()
            if r < 0.25 
                if random_bool() 
                    root = root + atom1 
                else
                    root = atom1 + root 
                end 
            elseif r < 0.5
                if random_bool()
                    root = root * atom2 
                else
                    root = atom2 * root
                end
            elseif r < 0.75
                if random_bool()
                    root = root + atom3 
                else
                    root = atom3 + root 
                end
            else
                if random_bool()
                    root = root + atom4 
                else
                    root = atom4 + root 
                end
            end
        end
    end

    root
end

function count_nodes_pre_ot(root::AbstractExpr)
    x = 0 

    UFL.@pre_order_traversal for s ∈ root 
        x += 1
    end

    x 
end

function count_nodes_post_ot(root::AbstractExpr)
    x = 0 

    UFL.@post_order_traversal for s ∈ root 
        x += 1
    end

    x
end

using BenchmarkTools

# Benchmarking Traversal:
# for n ∈ [100, 200, 300, 500, 1_000, 5_000, 10_000, 20_000] 
#     suite["$(n)_pre"] = @benchmarkable count_nodes_pre_ot(tree) setup=(tree=build_static_type_tree($n, $x, $u, $v))
#     suite["$(n)_post"] = @benchmarkable count_nodes_post_ot(tree) setup=(tree=build_static_type_tree($n, $x, $u, $v))
# end

function do_benchmarks()
    suite = BenchmarkGroup()
    
    mesh = UnitTriangleMesh()
    (x,y) = SpatialCoordinate(2, 2)
    
    V = VectorFunctionSpace(mesh, "CG", 1)
    
    u = UflFunction(V)
    v = TestFunction(V)
   
    for n ∈ [100, 200, 300, 500, 1_000, 5_000, 10_000, 20_000] 
        suite["$(n)_pre"] = @benchmarkable count_nodes_pre_ot(tree) setup=(tree=build_static_type_tree($n, $x, $u, $v))
        suite["$(n)_post"] = @benchmarkable count_nodes_post_ot(tree) setup=(tree=build_static_type_tree($n, $x, $u, $v))
    end

    tune!(suite)

    BenchmarkTools.run(suite, verbose=true, seconds=120)
end