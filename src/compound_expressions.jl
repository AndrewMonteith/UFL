export determinant_expr_2x2, determinant_expr_3x3, pseudo_determinant_expr

_det_2x2(A::AbstractExpr, i, j, k, l) = A[i, k] * A[j, l] - A[i, l] * A[j, k]
determinant_expr_2x2(A::AbstractExpr) = _det_2x2(A, 0, 1, 0, 1)

function codeterminant_expr_nxn(A::AbstractExpr, rows::Vector{Int}, cols::Vector{Int})::AbstractExpr
    length(rows) === 2 && _det_2x2(A, rows[1], rows[2], cols[1], cols[2])

    codet = 0.0
    r, subrows = rows[1], rows[2:end]

    for (i, c) ∈ enumerate(cols)
        subcols = [cols[i+1:end]..., cols[1:i-1]...]
        codet += A[r, c] * codeterminant_expr_nxn(A, subrows, subcols)
    end 

    codet
end

determinant_expr_3x3(A::AbstractExpr) = codeterminant_expr_nxn(A, [0, 1, 2], [0, 1, 2])