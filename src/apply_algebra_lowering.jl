export apply_algebra_lowering

lower_compound_algebra(expr::AbstractExpr, operands::VarTuple{AbstractExpr}) = reuse_if_untouched(expr, operands)

function lower_compound_algebra(::Determinant, operands::Tuple{AbstractExpr})::AbstractExpr
    expr = operands[1]
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

function lower_compound_algebra(::Transposed, op::Tuple{AbstractExpr})
    (i, j) = indices_n(2)

    as_tensor(op[1][i, j], (j, i))
end 

function lower_compound_algebra(::Trace, op::Tuple{AbstractExpr})::AbstractExpr
    i = Index() 
    op[1][i, i]
end

function lower_compound_algebra(::Dot, operands::Tuple{AbstractExpr, AbstractExpr})
    a, b = operands 

    ai = indices_n((length ∘ ufl_shape)(a)-1)
    bi = indices_n((length ∘ ufl_shape)(b)-1)
    k = Index()

    s = a[ai..., k] * b[k, bi...] 

    as_tensor(s, tuple(ai..., bi...))
end


apply_algebra_lowering(f::Union{Form, AbstractExpr}) = map_integrand_dags(lower_compound_algebra, f)