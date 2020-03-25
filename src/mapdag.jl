export map_expr_dag, map_integrand_dags, AbstractMapper, cached

abstract type AbstractMapper <: Function end

struct BaseMapper{V} 
    cached::Dict{AbstractExpr, V}
    BaseMapper() = new{AbstractExpr}(Dict{AbstractExpr, AbstractExpr}())
    BaseMapper{V}() where V = new{V}(Dict{AbstractExpr, V}())
end 

Base.getindex(m::AbstractMapper, x::AbstractExpr) = m.base.cached[x]
Base.getindex(m::AbstractMapper, x::Tuple{AbstractExpr, AbstractExpr})::Tuple{AbstractExpr, AbstractExpr} = (m[x[1]], m[x[2]])
Base.in(x::AbstractExpr, m::AbstractMapper) = x ∈ keys(m.base.cached) 

function (m::BaseMapper{AbstractExpr})(expr::AbstractExpr)::AbstractExpr 
    ops = ufl_operands(expr)
    operands = Tuple(m.cached[e] for e ∈ ops)
    ops === operands ? expr : reconstruct_expr(expr, operands)
end 

function map_expr_dag(mapper::AbstractMapper, expr::AbstractExpr)
    @UFL.post_order_traversal for node ∈ expr 
        node ∈ mapper && continue
        mapper.base.cached[node] = mapper(node)
    end

    mapper[expr]
end 

function map_integrands(func::F, form::Form) where F <: Function 
    mapped_integrals = [map_integrands(func, integral) for integral ∈ form.integrals]
    nonzero_integrals = filter(integral -> !(integral isa Zero), mapped_integrals)

    Form(tuple(nonzero_integrals...))
end

function map_integrands(func::F, integral::Integral)::Integral where F <: Function 
    reconstruct(integral; integrand=func(integral.integrand))
end

function map_integrands(func::F, expr::AbstractExpr)::AbstractExpr where F <: Function 
    func(expr)
end 

function map_integrand_dags(func::AbstractMapper, expr::Union{Form, Integral, AbstractExpr})
    map_integrands(e -> map_expr_dag(func, e), expr)
end 
