export apply_derivatives

abstract type DerivativeRulset <: Function end 

var_shape(g::DerivativeRulset) = g.base.var_shape
zero_terminal(g::DerivativeRulset, o::AbstractExpr) = Zero(tuple(ufl_shape(o)..., var_shape(g)...))

struct GenericDerivativeRuleset <: DerivativeRulset
    var_shape::DimensionTuple 
end 
var_shape(g::GenericDerivativeRuleset) = g.var_shape

(g::GenericDerivativeRuleset)(m::MultiIndexNode, _) = m
(g::GenericDerivativeRuleset)(c::AbstractConstantValue, _) = zero_terminal(g, c)
(g::GenericDerivativeRuleset)(s::Sum, ops::Tuple{AbstractExpr, AbstractExpr}) = ops[1] + ops[2]
function (g::GenericDerivativeRuleset)(p::Product, ops::Tuple{AbstractExpr, AbstractExpr})
    a, b = p.ufl_operands 
    da, db = ops 

    (da, db), ii = as_scalars(da, db)

    pa = Product(da, b)
    pb = Product(a, db)

    s = Sum(pa, pb)

    (!isempty(ii)) && (s = as_tensor(s, ii))
    
    s
end

function (g::GenericDerivativeRuleset)(p::Power, ops::Tuple{AbstractExpr, Zero})
    f, g = p.ufl_operands

    !is_true_scalar(f) && error("Expecting scalar f in f**g")
    !is_true_scalar(g) && error("Expecting scalar g in f**g")

    ops[1] * g * f^(g-1)
end 

function (g::GenericDerivativeRuleset)(p::Power, ops::Tuple{AbstractExpr, AbstractExpr})
    error("do not support generic power rule yet")
end

function (g::GenericDerivativeRuleset)(i::Indexed, ops::Tuple{AbstractExpr, AbstractExpr})
    Ap, ii = ops

    Ap === i.ufl_operands[1] && return i 

    if Ap isa ComponentTensor 
        # Unwrapping as_tensor(C[kk], jj)[ii] -> C[ll] 
        B, jj = ops[1].ufl_operands
        if B isa Indexed 
            C, kk = B.ufl_operands 

            if all(j ∈ kk for j ∈ jj)
                Cind = collect(jj.indices...)
                for (i, j) ∈ zip(ii.indices, jj.indices)
                    Cind[findfirst(x -> x === j, kk)] = i 
                end 

                return Indexed(C, MultiIndexNode(tuple(Cind...)))
            end 
        end
    end 


    r = (length ∘ ufl_shape)(Ap) - length(ii)
    
    if r > 0 
        kk = indices_n(r)
        as_tensor(Indexed(Ap, tuple(ii..., kk...)), kk)
    else
        Indexed(Ap, ii)
    end
end

struct GradRuleset <: DerivativeRulset
    base::GenericDerivativeRuleset
    id::Identity

    GradRuleset(gdim::Dimension) = new(GenericDerivativeRuleset((gdim,)), Identity(gdim))
end 

(g::GradRuleset)(sc::SpatialCoordinate, _) = g.id
(g::GradRuleset)(c::Constant, _) = is_cellwise_constant(c) ? zero_terminal(g) : Grad(c)
(g::GradRuleset)(m::AbstractExpr, x) = g.base(m, x)

derivative_dispatch(t::Terminal, op::Tuple{}) = t
derivative_dispatch(g::Grad, f::Tuple{AbstractExpr}) = map_expr_dag(f[1], GradRuleset(g.ufl_shape[end]))

derivative_dispatch(expr::AbstractExpr, operands::VarTuple{AbstractExpr}) =
    ufl_operands(expr) === operands ? expr : typeof(expr)(expr, operands)


apply_derivatives(expr::AbstractExpr) = map_expr_dag(expr, derivative_dispatch)