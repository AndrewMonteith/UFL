export map_expr_dag, remove_common_subexpressions, map_integrand_dags

function map_expr_dag(func::F, expr::AbstractExpr, ::Type{V}=Any) where {F<:Function, V<:Any}
    expr_vals = Dict{AbstractExpr, V}()
    vals = Dict{V, V}()
    
    @UFL.post_order_traversal for node ∈ expr 
        node ∈ keys(expr_vals) && continue

        v = func(node, map(op -> expr_vals[op], ufl_operands(node)))::V
        expr_vals[node] = get!(vals, v, v)
    end

    expr_vals[expr]
end 


function remove_common_subexpressions(root::AbstractExpr)
    seen_exprs = Dict{AbstractExpr, AbstractExpr}()

    remove_seen_expr(expr::Terminal, operands::Tuple{}) = expr 
    function remove_seen_expr(expr::Operator, operands::VarTuple{AbstractExpr})
        typeof(expr)(expr, operands)
    end 

    map_expr_dag(root, remove_seen_expr, AbstractExpr)
end


function map_integrands(func::F, form::Form) where F <: Function 
    mapped_integrals = [map_integrands(func, integral) for integral ∈ form.integrals]
    nonzero_integrals = filter(integral -> !(integral isa Zero), mapped_integrals)

    Form(tuple(nonzero_integrals...))
end

function map_integrands(func::F, integral::Integral) where F <: Function 
    reconstruct(integral; integrand=func(integral.integrand))
end

function map_integrands(func::F, expr::AbstractExpr)::AbstractExpr where F <: Function 
    func(expr)
end 

function map_integrand_dags(func::F, expr::Union{Form, Integral, AbstractExpr}) where F <: Function 
    map_integrands(e -> map_expr_dag(func, e), expr)
end 
