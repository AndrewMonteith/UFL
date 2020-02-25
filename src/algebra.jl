export Sum, Product, Divison, Power

function binary_show(io::IO, symbol::String, parent::Operator)
    ops = [parstr(parent, op) for op ∈ parent.ufl_operands]
    print(io, "$(ops[1]) $(symbol) $(ops[2])")
end 

@ufl_type struct Sum <: Operator 
    ufl_fields = (operands,)
    ufl_tags=(num_ops=2,)

    function Sum(a::AbstractExpr, b::AbstractExpr)
        shape = ufl_shape(a)
        if shape !== ufl_shape(b) 
            error("Cannot add expressions of different shape")
        end 

        free_indices = ufl_free_indices(a)
        if free_indices != ufl_free_indices(b) 
            error("Cannot add indices with different type")
        end 

        index_dimensions = ufl_index_dimensions(a)
        if index_dimensions != ufl_index_dimensions(b)
            error("Cannot add expressions with different index dimensions")
        end 

        if a isa Zero 
            return b 
        elseif b isa Zero 
            return a 
        end 

        a_is_scalar, b_is_scalar = a isa ScalarValue, b isa ScalarValue 

        if a_is_scalar && b_is_scalar
            return ScalarValue(a.val + b.val)
        elseif b_is_scalar
            a, b = b, a
        else
            # Do some kind of convoluted sorting algorithm which I hope I don't 
            # have to implement 
        end 

        new(@sig((a, b)))
    end
end

Base.show(io::IO, s::Sum) = binary_show(io, "+", s)
Base.:+(e1, e2) = Sum(as_ufl(e1), as_ufl(e2))
Base.:-(e1, e2) = Sum(as_ufl(e1), -as_ufl(e2))
Base.:-(e::AbstractExpr) = -1 * e

function merge_unqiue_indices(afi, afid, bfi, bfid)
    len_a, len_b = length(afi), length(bfi)

    len_a === 0 && return bfi, bfid 
    len_b === 0 && return afi, afid 

    ak, bk = 1, 1
    fi, fid = [], []

    while !(ak === len_a || bk === len_b)
        if afi[ak] < bfi[bk]
            push!(fi, afi[ak])
            push!(fid, afid[ak])
            ak += 1
        elseif afi[ak] > bfi[bk]
            push!(fi, bfi[bk])
            push!(fid, bfid[bk])
            bk += 1
        else
            push!(fi, afi[ak])
            push!(fid, afid[ak])
            ak += 1
            bk += 1
        end
    end

    if ak === len_a 
        if bk !== len_b 
            append!(fi, bfi[bk:end])
            append!(fid, bfid[bk:end])
        end
    elseif bk === len_b 
        if ak !== len_a 
            append!(fi, afi[ak:end])
            append!(fid, afid[ak:end])
        end
    end
    
    tuple(fi...), tuple(fid...) 
end

@ufl_type struct Product <: Operator 
    ufl_fields = (operands, free_indices, index_dimensions)
    ufl_tags=(nums_ops=2,)

    function Product(a::AbstractExpr, b::AbstractExpr) 
        (isempty ∘ ufl_shape)(a) || (isempty ∘ ufl_shape)(b) && error("product can only represent product of scalars")

        (a isa Zero || b isa Zero) && return Zero((), merge_unqiue_indices(ufl_free_indices(a), ufl_index_dimensions(a),
                                                                           ufl_free_indices(b), ufl_index_dimensions(b))...)

        a_is_scalar, b_is_scalar = a isa ScalarValue, b isa ScalarValue

        if a_is_scalar && b_is_scalar
            return ScalarValue(a.val * b.val)
        elseif a_is_scalar && a.val === 1
            return b 
        elseif b_is_scalar && b.val === 1
            return a 
        end 

        fi, fid = merge_unqiue_indices(ufl_free_indices(a), ufl_index_dimensions(a),
                                       ufl_free_indices(b), ufl_index_dimensions(b))

        new(@sig((a, b)), fi, fid)
    end
end
Base.show(io::IO, p::Product) = binary_show(io, "*", p)


function merge_overlappin_indices(afi, afid, bfi, bfid)
    len_a, len_b = length(afi), length(bfi)

    repeated_indices, repeated_index_dimensions = [], []
    
    for ii ∈ 1:len_a 
        for j ∈ 1:len_b 
            if afi[ii] == bfi[j]
                push!(repeated_indices, afi[ii])
                push!(repeated_index_dimensions, afid[ii])
                break
            end
        end
    end

    free_indices, index_dimensions = [], []

    for (i, d) ∈ zip((afi..., bfi...) , (afid..., bfid...))
        if i ∉ repeated_indices
            push!(free_indices, i)
            push!(index_dimensions, d)
        end
    end

    length(Set(free_indices)) === length(free_indices) || error("Free indices must not contain repeats")
    length(free_indices) + 2(length(repeated_indices)) === len_a+len_b || error("Expected only twice repeated indices")

    tuple(free_indices...), tuple(index_dimensions...), tuple(repeated_indices...), tuple(repeated_index_dimensions...)
end

function mult(a::AbstractExpr, b::AbstractExpr)
    fi, fid, ri, rid = merge_overlappin_indices(ufl_free_indices(a), ufl_index_dimensions(a),
                                                ufl_free_indices(b), ufl_index_dimensions(b))
   
    shape1, shape2 = ufl_shape(a), ufl_shape(b)
    rank1, rank2 = length(shape1), length(shape2)

    if rank1 === 0 && rank2 === 0 
        p = Product(a, b)
        ti = () 
    elseif rank1 === 0 || rank2 === 0 
        if rank2 === 0 
            a, b = b, a 
        end

        if a isa Zero || b isa Zero 
            shape = isempty(shape1) ? shape2 : shape1 
            return Zero(shape, fi, fid)
        end 

        ti = (indices_n ∘ length ∘ ufl_shape)(b)
        p = Product(a, b[ti...])
    elseif rank1 === 2 && (rank2 === 1 || rank2 === 2)
        !isempty(ri) && error("Not expecting repeate indices in non-scalar product.")

        (a isa Zero || b isa Zero) && return Zero(tuple(shape1[1:end-1]..., shape2[2:end]...), fi, fid)

        ai = indices_n(length(shape1) - 1)
        bi = indices_n(length(shape2) - 1)
        k = Index()

        p = a[ai..., k] * b[k, bi...]
        ti = tuple(ai..., bi...)
    else
        error("Invalid ranks $(rank1) and $(rank2) in product.")
    end

    if !isempty(ti)
        p = as_tensor(p, ti)
    end 

    for ii ∈ ri 
        p = IndexSum(p, (ii,))
    end 

    p
end
Base.:*(e1, e2) = mult(as_ufl(e1), as_ufl(e2))

is_true_scalar(a::AbstractExpr) = (isempty ∘ ufl_shape)(a) && (isempty ∘ ufl_free_indices)(a)

@ufl_type struct Division <: Operator 
    ufl_fields = (operands,)
    ufl_tags = (num_ops=2,)

    function Division(a::AbstractExpr, b::AbstractExpr)
        ufl_shape(a) !== () && error("expecting scalar numerator in Divison.")
        !is_true_scalar(b) && error("denominator must be a true scalar.")

        (a isa Zero || (b isa ScalarValue && b.val === 1)) && return a 

        (a isa ScalarValue && b isa ScalarValue) && ScalarValue(a.val / b.val)

        new(@sig((a, b)))
    end
end
Base.show(io::IO, div::Division) = binary_show(io, "/", div)

function _div(e1::AbstractExpr, e2::AbstractExpr)
    sh = ufl_shape(e1)
    
    if !isempty(sh) 
        ii = (indices_n ∘ length)(sh)
        d = Division(e1[ii...], e2)
        as_tensor(d, ii)
    else
        Divison(e1, e2)
    end
end
Base.:/(e1, e2) = _div(as_ufl(e1), as_ufl(e2))

@ufl_type struct Power <: Operator 
    ufl_fields = (operands,)
    ufl_tags = (num_ops=2,)

    function Power(a::AbstractExpr, b::AbstractExpr) 
        !is_true_scalar(a) && error("Cannot take the power of a non-scalar expression")
        !is_true_scalar(b) && error("Cnanot raise an expression to a non-scalar power")

        println("Creating a power of ", a, " ", b, " ", b isa ScalarValue)

        if a isa ScalarValue && b isa ScalarValue 
            return ScalarValue(a.val ^ b.val)
        elseif b isa Zero 
            return ScalarValue(1)
        elseif a isa Zero && b isa ScalarValue 
            b.val < 0 && error("Cannot raise 0 to a negative value")
            return Zero((), (), ())
        elseif b isa ScalarValue && b.val == 1
            return a
        end

        new(@sig((a, b)))
    end
end
Base.:^(e1, e2) = Power(as_ufl(e1), as_ufl(e2))
Base.show(io::IO, p::Power) = binary_show(io, "**", p)