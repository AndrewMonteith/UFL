export Sum, Product, Divison, Power

@ufl_type struct Sum <: Operator 
    ufl_fields = (operands,)
    ufl_tags = (inherit_shape_from_operand=1, inherit_indices_from_operand=1)
    
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

        new((a, b))
    end
end 

Base.:+(e1::AbstractExpr, e2::AbstractExpr) = Sum(e1, e2)
Base.:+(e1::AbstractExpr, e2::Real) = Sum(e1, as_ufl(e2))
Base.:+(e1::Real, e2::AbstractExpr) = Sum(as_ufl(e1), e2)
Base.:-(e1::AbstractExpr, e2::AbstractExpr) = Sum(e1, -e2)
Base.:-(e1::AbstractExpr, e2::Real) = Sum(e1, -as_ufl(e2))
Base.:-(e1::Real, e2::AbstractExpr) = Sum(as_ufl(e1), -e2)

Base.:-(e::AbstractExpr) = -1 * e

function merge_unqiue_indices(afi, afid, bfi, bfid)
    len_a, len_b = length(afi), length(bfi)

    len_a === 0 && return bfi, bfid 
    len_b === 0 && return afi, afid 

    ak, bk = 0, 0
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
    
    tuple(fi), tuple(fid) 
end

@ufl_type struct Product <: Operator 
    ufl_fields = (operands, free_indices, index_dimensions)

    function Product(a::AbstractExpr, b::AbstractExpr) 
        if ufl_shape(a) === () || ufl_shape(b) === () 
            error("product can only represent product of scalars")
        end

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

        new((a, b), merge_unqiue_indices(), merge_unqiue_indices(ufl_free_indices(a), ufl_index_dimensions(a),
                                                                 ufl_free_indices(b), ufl_index_dimensions(b))...)
    end
end

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

    length(Set(free_indices)) === length(free_indices) || error("Free indicies must not contain repeats")
    length(free_indices) + 2(length(repeated_indices)) === len_a+len_b || error("Expected only twice repeated indicies")

    tuple(free_indices), tuple(index_dimensions), tuple(repeated_indices), tuple(repeated_index_dimensions)
end

function mult(a::AbstractExpr, b::AbstractExpr)
    fi, fid, ri, rid = merge_overlappin_indices(ufl_free_indices(a), ufl_index_dimensions(a),
                                                ufl_free_indices(b), ufl_index_dimensions(b))
    
    shape1, shape2 = ufl_shape(a), ufl_shape(b)
    rank1, rank2 = length(shape1), length(shape2)

    if rank1 === 0 && rank2 === 0 
        p = Product(a, b)
    elseif rank1 === 0 || rank2 === 0 
        if rank2 === 0 
            a, b = b, a 
        end
        # rank2 === 0 && (a, b = b, a)

        if a isa Zero || b isa Zero 
            shape = shape1 === () ? shape2 : shape1 
            return Zero(shape, fi, fid)
        end 
    end

    error("not implemented")
end

is_true_scalar(a::AbstractExpr) = ufl_shape(a) === () && ufl_free_indices(a) === ()

@ufl_type struct Division <: Operator 
    ufl_fields = (operands,)
    ufl_tags = (inherit_indices_from_operand=0,)

    function Division(a::AbstractExpr, b::AbstractExpr)
        ufl_shape(a) !== () && error("expecting scalar numerator in Divison.")
        !is_true_scalar(a) && error("denominator must be a true scalar.")

        (a isa Zero || (b isa ScalarValue && b.val === 1)) && return a 

        (a isa ScalarValue && b isa ScalarValue) && ScalarValue(a.val / b.val)

        new((a, b))
    end
end
ufl_shape(d::Division) = ()

function _div(a::AbstractExpr, b::AbstractExpr)
    ufl_shape(a) !== () && error("cannot handle tensors yet")

    Divison(a, b)
end

Base.:/(e1::AbstractExpr, e2::AbstractExpr) = _div(e1, e2)
Base.:/(e1::AbstractExpr, e2::Real) = _div(e1, as_ufl(e2))
Base.:/(e1::Real, e2::AbstractExpr) = _div(as_ufl(e1), e2)

@ufl_type struct Power <: Operator 
    ufl_fields = (operands,)
    ufl_tags = (inherit_indices_from_operand=0,)

    function Power(a::AbstractExpr, b::AbstractExpr) 
        !is_true_scalar(a) && error("Cannot take the power of a non-scalar expression")
        !is_true_scalar(b) && error("Cnanot raise an expression to a non-scalar power")

        if a isa ScalarValue && b isa ScalarValue 
            ScalarValue(a.val ^ b.val)
        elseif b isa Zero 
            ScalarValue(1)
        elseif a isa Zero && b isa ScalarValue 
            b.val < 0 && error("Cannot raise 0 to a negative value")
            Zero((), (), ())
        elseif b isa ScalarValue && b.val === 1
            a
        else
            new((a, b))
        end
    end
end
ufl_shape(p::Power) = ()

