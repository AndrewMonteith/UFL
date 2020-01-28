export ComponentTensor, as_tensor

function remove_indices(fi::MultiIndex, fid::DimensionTuple, rfi::MultiIndex)
    isempty(rfi) && return fi, fid 

    rfip = sort(zip(ind, i) for (i, ind) ∈ enumerate(rfi))
    rfi_len, fi_len = length(rfi), length(fi) 

    shape = collect(-1 for _ ∈ 1:rfi_len)
    newfiid::Array{Tuple{AbstractIndex, Dimension}} = []

    k, pos = 1, 1 

    while pos < fi_len
        rk = rfip[k][1]

        while fi[pos] < rk 
            push!(newfiid, (fi[pos], fid[pos]))
            pos += 1
        end 

        removed = 0
        while pos < nfi && fi[pos] === rk 
            shape[rfip[k][2]] = fid[pos]
            pos += 1
            removed += 1
        end 

        # We expect to find each index from rfi in fi
        removed === 0 && error("Index to be removed $(rk) not part of indices $(fi)")

        k += 1
        if k > rfi_len 
            pos <= fi_len && append!(newfiid, zip(fi[pos:end-1], fid[pos:end-1]))
            break 
        end 
    end

    fi, fid = isempty(newfiid) ? ((), ()) : zip(newfiid...)
    
    fi, fid, tuple(shape)
end

@ufl_type struct ComponentTensor <: Operator 
    ufl_fields = (operands,)

    function ComponentTensor(expr::AbstractExpr, indices::VarTuple{Index})
        ufl_shape(expr) !== () && error("Expecting scalar valued expression.")

        fi, fid, sh = remove_indices(ufl_free_indices(expr),
                                     ufl_index_dimensions(expr),
                                     indices)

        new(sh, fi, fid, (expr, indices))
    end 
end 

@ufl_type struct ListTensor <: Operator 
    ufl_fields = (operands,)

    function ListTensor(exprs::Union{VarTuple{AbstractExpr}, AbstractArray{AbstractExpr}})
        e_shape = ufl_shape(exprs[1])
        e_fi = ufl_free_indices(exprs[1])
        e_fid = ufl_index_dimensions(exprs[1])

        all(e -> ufl_shape(e) == e_shape, exprs) || error("All subexpressions must have the same shape")
        all(e -> ufl_free_indices(e) == e_fi) || error("All components must have same free indices")
        all(e -> ufl_index_dimensions(e) == e_fid) || error("All component shave different free index dimensions")

        exprs = exprs isa Tuple ? exprs : tuple(exprs...)

        new((), (), (), exprs)
    end
end

ufl_shape(lt::ListTensor) = tuple(length(lt.ufl_operands), ufl_shape(lt.ufl_operands[1])...)

as_tensor(expr::AbstractExpr) = expr 
function as_tensor(exprs::Union{VarTuple{AbstractExpr}, AbstractArray{AbstractExpr}})
    recursive_convert(expr::Union{Tuple, AbstractArray}) = ListTensor(collect(recursive_convert(e) for e ∈ expr)...)
    recursive_convert(expr) = as_ufl(expr)

    recursive_convert(exprs)
end

function as_tensor(expr::AbstractExpr, indices::MultiIndex)
    isempty(indices) && return expr 

    if expr isa Indexed 
        A, ii = ufl_operands(expr)
        indices == ii && return A 
    end 

    ComponentTensor(expr, indices)
end
as_tensor(expr::AbstractExpr, index::AbstractExpr) = as_tensor(expr, (index,))

function as_matrix(expr::AbstractExpr)
    length(ufl_shape(expr)) !== 2 && error("Expecting rank 2 tensor.")
    expr 
end

as_matrix(exprs::Union{VarTuple{AbstractExpr}, AbstractArray{AbstractExpr}}) = as_tensor(exprs)
function as_matrix(exprs::Union{VarTuple{AbstractExpr}, AbstractArray{AbstractExpr}}, indices::MultiIndex) 
    length(indices) !== 2 && error("Expecting exactly two indices")

    as_tensor(exprs, indices)
end