export apply_derivatives, reuse_if_untouched

abstract type DerivativeRulset <: Function end 

reuse_if_untouched(expr::AbstractExpr, operands::VarTuple{AbstractExpr})::AbstractExpr = 
    ufl_operands(expr) === operands ? expr : typeof(expr)(expr, operands)

var_shape(g::DerivativeRulset) = g.base.var_shape
zero_terminal(g::DerivativeRulset, o::AbstractExpr) = Zero(tuple(ufl_shape(o)..., var_shape(g)...))
zero_operator(g::DerivativeRulset, o::AbstractExpr) = Zero(tuple(ufl_shape(o)..., var_shape(g)...), ufl_free_indices(o), ufl_index_dimensions(o))

struct GenericDerivativeRuleset <: DerivativeRulset
    var_shape::DimensionTuple 
end 
var_shape(g::GenericDerivativeRuleset) = g.var_shape

(g::GenericDerivativeRuleset)(m::MultiIndexNode, _) = m
(g::GenericDerivativeRuleset)(c::AbstractConstantValue, _) = zero_terminal(g, c)
(g::GenericDerivativeRuleset)(s::Sum, ops::Tuple{AbstractExpr, AbstractExpr}) = ops[1] + ops[2]
(g::GenericDerivativeRuleset)(::IndexSum, ops::Tuple{AbstractExpr, AbstractExpr}) = IndexSum(ops...)
(g::GenericDerivativeRuleset)(f::UflFunction, ops::Tuple{}) = reuse_if_untouched(f, ops)

function (g::GenericDerivativeRuleset)(ln::Ln, op::Tuple{AbstractExpr})
    f, = ufl_operands(ln) 
    
    f isa Zero && error("Division by zero")
    
    op[1]/f
end 

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
    elseif r === 0
        Indexed(Ap, ii)
    else
        println(i)
        println(ops)
        println("Expression:", Ap)
        println("ii:", ii)
        error("wat")
    end
end

(g::GenericDerivativeRuleset)(ct::CompoundTensorOperator, ops::Tuple{Zero, AbstractExpr}) = zero_operator(ct)
function (g::GenericDerivativeRuleset)(ct::ComponentTensor, ops::Tuple{AbstractExpr, AbstractExpr})
    Ap, ii = ops 
    Ap, jj = as_scalars(Ap) # Ap is tuple because i've only implemented group as_scalars

    as_tensor(Ap[1], tuple(ii.indices..., jj...))
end



struct GradRuleset <: DerivativeRulset
    base::GenericDerivativeRuleset
    id::Identity

    GradRuleset(gdim::Dimension) = new(GenericDerivativeRuleset((gdim,)), Identity(gdim))
end 

(g::GradRuleset)(sc::SpatialCoordinate, _) = g.id
(g::GradRuleset)(c::Union{Constant, UflFunction}, _) = is_cellwise_constant(c) ? zero_terminal(g) : Grad(c)
(g::GradRuleset)(m::AbstractExpr, x) = g.base(m, x)


struct GateauxDerivativeRuleset <: DerivativeRulset
    base::GenericDerivativeRuleset
    w::VarTuple{AbstractExpr} # coefficients 
    v::VarTuple{AbstractExpr} # arguments 
    w2v::Dict{AbstractExpr, AbstractExpr} # coefficients -> arguments
    cd::Dict{AbstractExpr, AbstractExpr} # nonzero df/dw

    function GateauxDerivativeRuleset(coefficients::ExprList, arguments::ExprList, coefficient_derivatives::ExprList)
        w2v = Dict(w => v for (w, v) ∈ zip(coefficients, arguments))

        cd = ufl_operands(coefficient_derivatives)
        cd′ = if isempty(cd)
            Dict{AbstractExpr, AbstractExpr}()
        else
            Dict(cd[2*i] => cd[2*i+1] for i ∈ 1:(length(cd)÷2))
        end      

        new(GenericDerivativeRuleset(()), ufl_operands(coefficients), ufl_operands(arguments), w2v, cd′)
    end 
end 
function (g::GateauxDerivativeRuleset)(o::Coefficient, ops::Tuple{})
    _do = get(g.w2v, o, nothing)

    _do !== nothing && return _do 

    dos = get(g.cd, o, nothing)
    if dos === nothing 
        Zero(ufl_shape(o))
    else
        error("what are we doing here")
    end
end 

(g::GateauxDerivativeRuleset)(o::Argument, ops::Tuple{}) = zero_terminal(g, o)
(g::GateauxDerivativeRuleset)(x::AbstractExpr, ops::VarTuple{AbstractExpr}) = g.base(x, ops)
function (g::GateauxDerivativeRuleset)(grad::Grad, _::VarTuple{AbstractExpr})     

    # Need to handle: grad(tr(grad(u)))

    # Count how deep the grad is
    ngrads, o = 0, grad
    while o isa Grad 
        o, = ufl_operands(o)
        ngrads += 1
    end 

    o isa AbstractFormArgument || error("grad must contains form argument")

    for (w, v) ∈ zip(g.w, g.v)
        if o == w && v isa AbstractFormArgument
            for i in 1:ngrads 
                v = Grad(v)
            end 

            return v 
        end 
    end 

    error("not implemented this yet")
end


derivative_dispatch(t::Terminal, op::Tuple{}) = t
derivative_dispatch(g::Grad, f::Tuple{AbstractExpr}) = map_expr_dag(GradRuleset(g.ufl_shape[end]), f[1])
derivative_dispatch(expr::AbstractExpr, operands::VarTuple{AbstractExpr}) = reuse_if_untouched(expr, operands)
function derivative_dispatch(cd::CoefficientDerivative, ops::NTuple{4, AbstractExpr})
    dummpy, w, v, cd = ufl_operands(cd) 
    rules = GateauxDerivativeRuleset(w, v, cd)
    map_expr_dag(rules, ops[1])
end 

apply_derivatives(expr::Form) = map_integrand_dags(derivative_dispatch, expr)
apply_derivatives(expr::AbstractExpr) = map_expr_dag(derivative_dispatch, expr)