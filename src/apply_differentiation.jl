export apply_derivatives, reuse_if_untouched

var_shape(g::AbstractMapper) = g.var_shape
zero_terminal(g::AbstractMapper, o::AbstractExpr) = Zero(tuple(ufl_shape(o)..., var_shape(g)...))
zero_operator(g::AbstractMapper, o::AbstractExpr) = Zero(tuple(ufl_shape(o)..., var_shape(g)...), ufl_free_indices(o), ufl_index_dimensions(o))

function indexed_derivative(m::AbstractMapper, i::Indexed)::AbstractExpr
    Ap, ii = m[ufl_operands(i)]

    Ap === i.ufl_operands[1] && return i 

    if Ap isa ComponentTensor 
        # Unwrapping as_tensor(C[kk], jj)[ii] -> C[ll] 
        B, jj = ufl_operands(Ap)
        if B isa Indexed 
            C, kk = ufl_operands(B)

            if all(j ∈ kk for j ∈ jj)
                Cind = collect(kk.indices)
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

generic_derivative_rule(mapper::AbstractMapper, m::MultiIndexNode) = m 
generic_derivative_rule(mapper::AbstractMapper, c::AbstractConstantValue) = zero_terminal(mapper, c)
generic_derivative_rule(mapper::AbstractMapper, f::UflFunction) = mapper.base(f)
generic_derivative_rule(mapper::AbstractMapper, i::Indexed) = indexed_derivative(mapper, i)
generic_derivative_rule(mapper::AbstractMapper, c::Conj) = conj(mapper[c.ufl_operands[1]])

function generic_derivative_rule(mapper::AbstractMapper, s::Sum)
    a, b = mapper[ufl_operands(s)]
    a + b
end 

function generic_derivative_rule(mapper::AbstractMapper, i::IndexSum)
    a, b = mapper[ufl_operands(i)]
    IndexSum(a, b)
end

function generic_derivative_rule(mapper::AbstractMapper, l::Ln)
    f = ufl_operands(l)[1] 
    f′ = mapper[f]

    f′  isa Zero && error("Division by zero")

    f′/f 
end 

function generic_derivative_rule(mapper::AbstractMapper, p::Product)
    a, b = ufl_operands(p)
    a′, b′ = mapper[a], mapper[b]

    (da, db), ii = as_scalars(a′, b′)
    
    pa = Product(da, b)
    pb = Product(a, db)
    
    s = Sum(pa, pb)
    
    (!isempty(ii)) && (s = as_tensor(s, ii))
    
    s
end

function generic_derivative_rule(mapper::AbstractMapper, p::Power)
    f, g = ufl_operands(p)

    !is_true_scalar(f) && error("Expecting scalar f in f**g")
    !is_true_scalar(g) && error("Expecting scalar g in f**g")

    f′, g′ = mapper[f], mapper[g]

    if g′ isa Zero
        f′ * g * f^(g-1)
    else
        error("do not support generic power rule yet")
    end
end

function generic_derivative_rule(mapper::AbstractMapper, ct::ComponentTensor)
    Ap, ii = mapper[ufl_operands(ct)]

    if Ap isa Zero 
        zero_operator(mapper, ct)
    else
        Ap, jj = as_scalars(Ap) # Ap is tuple because i've only implemented group as_scalars
        as_tensor(Ap[1], tuple(ii.indices..., jj...))
    end
end

function generic_derivative_rule(mapper::AbstractMapper, d::Division)
    f, g = ufl_operands(d)
    f′, g′ = mapper[f], mapper[g]

    !UFL.is_scalar(f) && error("Not expecting nonscalar numerator")
    !UFL.is_scalar(g) && error("Not expecting nonscalar denominator")
    
    sd, si = as_scalar(d)
    sgp, gi = as_scalar(g′)
    
    o_gp = sd*sgp
    if !(isempty(si) && isempty(gi))
        o_gp = as_tensor(o_gp, tuple(si..., gi...))
    end 
    
    (f′-o_gp)/g
end


struct GradDerivativeMapper <: AbstractMapper 
    base::BaseMapper
    var_shape::DimensionTuple
    id::Identity

    GradDerivativeMapper(gdim::Dimension) = new(BaseMapper(), (gdim,), Identity(gdim))
end

(g::GradDerivativeMapper)(sc::SpatialCoordinate) = g.id
(g::GradDerivativeMapper)(c::Union{Constant, UflFunction}) = is_cellwise_constant(c) ? zero_terminal(g, c) : Grad(c)
(g::GradDerivativeMapper)(arg::Argument) = Grad(arg)
function (g::GradDerivativeMapper)(gr::Grad)
    (typeof(gr.ufl_operands[1]) isa Grad || typeof(gr.ufl_operands[1]) <: Terminal) || error("Expecting only grads applied to a terminal")
    Grad(gr)
end
(g::GradDerivativeMapper)(x::AbstractExpr) = generic_derivative_rule(g, x)


struct GateauxDerivativeMapper <: AbstractMapper 
    base::BaseMapper
    var_shape::DimensionTuple
    w::VarTuple{AbstractExpr} # coefficients 
    v::VarTuple{AbstractExpr} # arguments 
    w2v::Dict{AbstractExpr, AbstractExpr} # coefficients -> arguments
    cd::Dict{AbstractExpr, AbstractExpr} # nonzero df/dw

    function GateauxDerivativeMapper(coefficients::ExprList, arguments::ExprList, coefficient_derivatives::ExprList)
        w2v = Dict(w => v for (w, v) ∈ zip(coefficients, arguments))

        cd = ufl_operands(coefficient_derivatives)
        cd′ = if isempty(cd)
            Dict{AbstractExpr, AbstractExpr}()
        else
            Dict(cd[2*i] => cd[2*i+1] for i ∈ 1:(length(cd)÷2))
        end      

        new(BaseMapper(), (), ufl_operands(coefficients), ufl_operands(arguments), w2v, cd′)
    end 
end 

(g::GateauxDerivativeMapper)(o::Argument) = zero_terminal(g, o)
(g::GateauxDerivativeMapper)(x::AbstractExpr) = generic_derivative_rule(g, x)

function (g::GateauxDerivativeMapper)(o::Coefficient)::AbstractExpr 
    _do = get(g.w2v, o, nothing)

    _do !== nothing && return _do 

    dos = get(g.cd, o, nothing)
    if dos === nothing 
        Zero(ufl_shape(o))
    else
        error("what are we doing here")
    end
end

function (g::GateauxDerivativeMapper)(grad::Grad)::AbstractExpr 
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

struct DerivativeMapper <: AbstractMapper 
    base::BaseMapper 
    DerivativeMapper() = new(BaseMapper())
end 

(d::DerivativeMapper)(t::Terminal) = t 
(d::DerivativeMapper)(g::Grad) = map_expr_dag(GradDerivativeMapper(g.ufl_shape[end]), d[ufl_operands(g)[1]])
(d::DerivativeMapper)(i::Indexed) = indexed_derivative(d, i)
(d::DerivativeMapper)(expr::AbstractExpr) = d.base(expr)

function (d::DerivativeMapper)(cd::CoefficientDerivative)
    dummy, w, v, cd = ufl_operands(cd) 
    rules = GateauxDerivativeMapper(w, v, cd)
    map_expr_dag(rules, d[dummy])
end

apply_derivatives(expr::Union{Form, AbstractExpr}) = map_integrand_dags(DerivativeMapper(), expr)