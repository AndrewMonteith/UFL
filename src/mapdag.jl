export map_expr_dag

function map_expr_dag(expr::AbstractExpr, func::F, ::Type{V}=Any) where {F<:Function, V<:Any}
    expr_vals = Dict{AbstractExpr, V}()
    vals = Dict{V, V}()
    
    @UFL.post_order_traversal for node ∈ expr 
        node ∈ keys(expr_vals) && continue

        v = func(node, map(op -> expr_vals[op], ufl_operands(node)))::V
        expr_vals[node] = get!(vals, v, v)
    end

    expr_vals[expr]
end 

function subexpression_replacement(root::AbstractExpr)
    seen_exprs = Dict{AbstractExpr, AbstractExpr}()

    function map_dag(expr::Terminal, operands::Tuple{})
        get!(seen_exprs, expr, expr)
    end
end