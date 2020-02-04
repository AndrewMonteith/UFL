export ComponentTensor, as_tensor, as_matrix

function remove_indices(fi::MultiIndex, fid::DimensionTuple, rfi::MultiIndex)
    isempty(rfi) && return fi, fid 

    rfip = (sort ∘ collect)((ind, i) for (i, ind) ∈ enumerate(rfi))
    rfi_len, fi_len = length(rfi), length(fi) 

    shape = collect(-1 for _ ∈ 1:rfi_len)
    newfiid::Array{Tuple{AbstractIndex, Dimension}} = []

    k, pos = 1, 1 

    while pos <= fi_len
        rk = rfip[k][1]

        while fi[pos] < rk 
            push!(newfiid, (fi[pos], fid[pos]))
            pos += 1
        end 

        removed = 0
        while pos <= fi_len && fi[pos] === rk 
            shape[rfip[k][2]] = fid[pos]
            pos += 1
            removed += 1
        end 

        # We expect to find each index from rfi in fi
        removed === 0 && error("Index to be removed $(rk) not part of indices $(fi)")

        k += 1
        if k > rfi_len 
            pos <= fi_len && append!(newfiid, zip(fi[pos:end], fid[pos:end]))
            break 
        end 
    end

    fi, fid = isempty(newfiid) ? ((), ()) : zip(newfiid...)
    
    fi, fid, tuple(shape...)
end

@ufl_type struct ComponentTensor <: Operator 
    ufl_fields = (operands,)

    function ComponentTensor(expr::AbstractExpr, indices::VarTuple{Index})
        ufl_shape(expr) !== () && error("Expecting scalar valued expression.")

        fi, fid, sh = remove_indices(ufl_free_indices(expr),
                                     ufl_index_dimensions(expr),
                                     indices)

        new(sh, fi, fid, (expr, as_ufl(indices)))
    end 
end 
Base.show(io::IO, ct::ComponentTensor) = print(io, "{ A | A_{$(ct.ufl_operands[2]) = $(ct.ufl_operands[1])}")

@ufl_type struct ListTensor <: Operator 
    ufl_fields = (operands,)

    function ListTensor(exprs::Union{VarTuple, AbstractArray})
        e_shape = ufl_shape(exprs[1])
        e_fi = ufl_free_indices(exprs[1])
        e_fid = ufl_index_dimensions(exprs[1])

        all(e -> ufl_shape(e) == e_shape, exprs) || error("All subexpressions must have the same shape")
        all(e -> ufl_free_indices(e) == e_fi, exprs) || error("All components must have same free indices")
        all(e -> ufl_index_dimensions(e) == e_fid, exprs) || error("All component shave different free index dimensions")

        exprs = exprs isa Tuple ? exprs : tuple(exprs...)

        new((), (), (), exprs)
    end
end

function Base.getindex(lt::ListTensor, key...)
    key_i = indices(key)

    index(d::Dimension) = true, d 
    index(i::FixedIndex) = true, i.d 
    index(i) = false, -1 

    valid, k = index(key_i[1])
    if valid 
        op = ufl_operands(lt)[k]

        return length(key_i) === 1 ? op : op[key[2:end]...]
    else
        invoke(Base.getindex, Tuple{AbstractExpr, Vararg}, lt, key...)
    end
end

ufl_shape(lt::ListTensor) = tuple(length(lt.ufl_operands), ufl_shape(lt.ufl_operands[1])...)

as_tensor(expr::AbstractExpr) = expr 
function as_tensor(exprs::Union{VarTuple, AbstractArray})
    recursive_convert(expr::Union{Tuple, AbstractArray}) = ListTensor(collect(recursive_convert(e) for e ∈ expr))
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