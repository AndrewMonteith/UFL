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

Base.show(io::IO, t::Transposed) = "$(parstr(t, (first ∘ ufl_operands)(t)))^T"

@ufl_type struct Trace <: CompoundTensorOperator 
    ufl_fields = (operands,)
    ufl_tags = (num_ops=1,)

    function Trace(A::AbstractExpr)
        (length ∘ ufl_shape)(A) !== 2 && error("Trace of tensor with rank != 2 is undefined")
        new(@sig((A,)))
    end 

    Trace(z::Zero) = Zero((), ufl_free_indices(z), ufl_index_dimensions(z))
end 

Base.show(io::IO, t::Trace) = "tr($((first ∘ ufl_operands)(t)))"


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