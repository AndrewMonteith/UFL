export tr, det

abstract type CompoundTensorOperator <: Operator end 

@ufl_type struct Transposed <: CompoundTensorOperator
    ufl_fields = (operands,shape)
    ufl_tags = (num_ops=1,)

    function Transposed(z::Zero)
        a, b = ufl_shape(z)
        Zero((b,a), ufl_free_indices(z), ufl_index_dimensions(z))
    end

    function Transposed(x::AbstractExpr)
        sh = ufl_shape(x)
        length(sh) !== 2 && error("Transposed only defined for rank 2 tensors")
        a, b = sh 
        new(@sig((x,)), (b, a))
    end
end 

Base.show(io::IO, t::Transposed) = print(io, "$(parstr(t, (first ∘ ufl_operands)(t)))^T")

@ufl_type struct Trace <: CompoundTensorOperator 
    ufl_fields = (operands,)
    ufl_tags = (num_ops=1,)

    function Trace(A::AbstractExpr)
        (length ∘ ufl_shape)(A) !== 2 && error("Trace of tensor with rank != 2 is undefined")
        new(@sig((A,)))
    end 

    Trace(z::Zero) = Zero((), ufl_free_indices(z), ufl_index_dimensions(z))
end 

Base.show(io::IO, t::Trace) = print(io, "tr($((first ∘ ufl_operands)(t)))")
tr(x::AbstractExpr) = Trace(x)


@ufl_type struct Determinant <: CompoundTensorOperator 
    ufl_fields = (operands,shape)
    ufl_tags=(num_ops=1,)

    Determinant(z::Zero) = Zero((), ufl_free_indices(z), ufl_index_dimensions(z))

    function Determinant(A::AbstractExpr)
        sh = ufl_shape(A)
        r = length(sh)
        Afi = ufl_free_indices(A)

        (r !== 0 && r !== 2) && error("Determinant of tensor with rank != 2 is undefined")
        (r === 2 & sh[1] !== sh[2]) && error("Cannot take determinant of rectangular rank 2 tensor")
        (!isempty(Afi)) && error("Not expecting free indices in determinant")

        r === 0 && return A 

        new(@sig((A,)), sh)
    end
end

Base.show(io::IO, d::Determinant) = print(io, "det($((first ∘ ufl_operands)(d)))")
det(x::AbstractExpr) = Determinant(x)


# @inline _getproperty(x::AbstractExpr, ::Val{s}) where s = getfield(x, s)
# @inline _getproperty(x::AbstractExpr, ::Val{:T}) = Transposed(x)

# @inline Base.getproperty(x::AbstractExpr, s::Symbol) = _getproperty(x, Val{s}())
# Intercepted this allows us to have transposed operators but is a performance bottleneck.
# Differentiation 6.2ms -> 8.4ms 
# Building Tree 404 microseconds -> 3ms 
function Base.getproperty(x::AbstractExpr, s::Symbol)
    if s === :T 
        Transposed(x)
    else
        getfield(x, s)
    end 
end 