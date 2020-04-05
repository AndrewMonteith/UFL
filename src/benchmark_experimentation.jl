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

function build_diff_example_1(x::AbstractExpr, u::AbstractExpr, v::AbstractExpr)
    root = x
    for i ∈ 2:10
        if random_bool()
            root = root + x^i 
        else
            root = x^i + root 
        end 
    end
    root
end

function build_diff_example_2(x::AbstractExpr, u::AbstractExpr, v::AbstractExpr)
    root = x - 1/x
    for i ∈ 2:10 
        term = (i*x^i - i/x^i)
        if random_bool()
            root = root + term
        else
            root = term + root 
        end 
    end
    root
end

function build_diff_example_3(x::AbstractExpr, u::AbstractExpr, v::AbstractExpr)
    root = 0
    for i ∈ 1:10 
        term = x^i*ln(x) + 1/(x^i*ln(x))
        if random_bool()
            root = root + term
        else
            root = term + root 
        end 
    end
    root
end

function build_diff_example_4(x::AbstractExpr, u::AbstractExpr, v::AbstractExpr)
    root = 0
    for i ∈ 1:10 
        term = x^i*dot(u, v) + ln(1/x^i)*inner(u, v)
        if random_bool()
            root = root + term
        else
            root = term + root 
        end 
    end
    root
end

function build_diff_example_5(x::AbstractExpr, u::AbstractExpr, v::AbstractExpr)
    root = 0
    for i ∈ 1:10 
        term = x^i*det(grad(u)) + ln(1/x^i)*tr(grad(u))
        if random_bool()
            root = root + term
        else
            root = term + root 
        end 
    end
    root
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
    
    V = VectorFunctionSpace(mesh, "CG", 1);
    
    u = UflFunction(V)
    v = TestFunction(V)
  
    suite["example_1"] = @benchmarkable apply_derivatives(grad(tree)) setup=(tree=UFL.apply_algebra_lowering(build_diff_example_1($x, $u, $v)))
    suite["example_2"] = @benchmarkable apply_derivatives(grad(tree)) setup=(tree=UFL.apply_algebra_lowering(build_diff_example_2($x, $u, $v)))
    suite["example_3"] = @benchmarkable apply_derivatives(grad(tree)) setup=(tree=UFL.apply_algebra_lowering(build_diff_example_3($x, $u, $v)))
    suite["example_4"] = @benchmarkable apply_derivatives(grad(tree)) setup=(tree=UFL.apply_algebra_lowering(build_diff_example_4($x, $u, $v)))
    suite["example_5"] = @benchmarkable apply_derivatives(grad(tree)) setup=(tree=UFL.apply_algebra_lowering(build_diff_example_5($x, $u, $v)))
    # for n ∈ [100, 200, 300, 500, 1_000, 5_000, 10_000, 20_000] 
    #     suite["$(n)"] = @benchmarkable UFL.apply_derivatives(grad(tree)) setup=(tree=UFL.apply_algebra_lowering(build_static_type_tree($n, $x, $u, $v)))
    # end

    tune!(suite)

    BenchmarkTools.run(suite, verbose=true, seconds=8)
end