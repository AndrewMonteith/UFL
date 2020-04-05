export apply_algebra_lowering

struct LowerAlgebra <: AbstractMapper 
    base::BaseMapper

    LowerAlgebra() = new(BaseMapper())
end
(l::LowerAlgebra)(x::AbstractExpr) = l.base(x)

function (l::LowerAlgebra)(d::Determinant)
    o = ufl_operands(d)

    expr = l[o[1]]
    sh = ufl_shape(expr)

    if sh === () 
        expr 
    elseif sh[1] === sh[2]
        if sh[1] === 1
            expr[0, 0]
        elseif sh[1] === 2
            determinant_expr_2x2(expr)
        elseif sh[1] === 3
            determinant_expr_3x3(expr)
        end
    else
        error("we don't support this for now")
    end
end

function (l::LowerAlgebra)(t::Transposed)
    (i, j) = indices_n(2)
    op = l[ufl_operands(t)[1]] # cached(l, ufl_operands(t)[1])

    as_tensor(op[i, j], (j, i))
end 

function (l::LowerAlgebra)(t::Trace)::AbstractExpr
    i = Index() 
    l[ufl_operands(t)[1]][i, i]
end

function (l::LowerAlgebra)(d::Dot)
    a, b = l[ufl_operands(d)]

    ai = indices_n((length ∘ ufl_shape)(a)-1)
    bi = indices_n((length ∘ ufl_shape)(b)-1)
    k = Index()

    s = a[ai..., k] * b[k, bi...] 

    as_tensor(s, tuple(ai..., bi...))
end

function (l::LowerAlgebra)(i::Inner)
    a, b = l[ufl_operands(i)]

    ash, bsh = ufl_shape(a), ufl_shape(b)
    
    ash !== bsh && error("Mismatching shape") 

    ii = (indices_n ∘ length)(ash)

    return a[ii...] * conj(b[ii...])
end

apply_algebra_lowering(f::Union{Form, AbstractExpr}) = map_integrand_dags(LowerAlgebra(), f)




# lower_compound_algebra(expr::AbstractExpr, cached::Dict{AbstractExpr, AbstractExpr}) = reuse_if_untouched_(expr, cached)

# function lower_compound_algebra(d::Determinant, cached::Dict{AbstractExpr, AbstractExpr})::AbstractExpr
#     o = ufl_operands(d)

#     expr = cached[o[1]]
#     sh = ufl_shape(expr)

#     if sh === () 
#         expr 
#     elseif sh[1] === sh[2]
#         if sh[1] === 1
#             expr[0, 0]
#         elseif sh[1] === 2
#             determinant_expr_2x2(expr)
#         elseif sh[1] === 3
#             determinant_expr_3x3(expr)
#         end
#     else
#         error("we don't support this for now")
#     end
# end

# function lower_compound_algebra(t::Transposed, cached::Dict{AbstractExpr, AbstractExpr})
#     (i, j) = indices_n(2)
#     op = cached[ufl_operands(t)[1]]

#     as_tensor(op[i, j], (j, i))
# end 

# function lower_compound_algebra(t::Trace, cached::Dict{AbstractExpr, AbstractExpr})::AbstractExpr
#     i = Index() 
#     cached[ufl_operands(t)[1]][i, i]
# end

# function lower_compound_algebra(d::Dot, cached::Dict{AbstractExpr, AbstractExpr})
#     x, y = ufl_operands(d) 
#     a, b = cached[x], cached[y]

#     ai = indices_n((length ∘ ufl_shape)(a)-1)
#     bi = indices_n((length ∘ ufl_shape)(b)-1)
#     k = Index()

#     s = a[ai..., k] * b[k, bi...] 

#     as_tensor(s, tuple(ai..., bi...))
# end


# apply_algebra_lowering(f::Union{Form, AbstractExpr}) = map_integrand_dags(lower_compound_algebra, f)