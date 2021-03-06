export Zero, Identity, ScalarValue, as_ufl

abstract type AbstractConstantValue <: Terminal end 

is_cellwise_constant(::AbstractConstantValue) = true
ufl_domains(::AbstractConstantValue) = ()

@ufl_type struct Zero <: AbstractConstantValue
    ufl_fields = (shape, free_indices, index_dimensions)

    function Zero(shape::DimensionTuple=())
        new(@sig(shape), @sig(()), @sig(()))
    end

    """
        Constructor for the old format 
        Example old format code:
            i = Index(2)
            j = Index(3)

            z = Zero((), (j, i). Dict{i => 3, j => 5})
    """
    function Zero(shape::DimensionTuple, free_indices::VarTuple{Index}, index_dimensions::Dict{Index, Dimension})
        sorted_indices = sort(collect(free_indices); by=x -> x.id)
        dimensions = Tuple(index_dimensions[index] for index in sorted_indices)
        new(@sig(shape), @sig(Tuple(sorted_indices)), @sig(index_dimensions))
    end

    """
        Constructor for the new format 
        Example new format code:
            z = Zero((), (2, 4), (3, 5))
    """
    function Zero(shape::DimensionTuple, free_indices::VarTuple{Index}, index_dimensions::VarTuple{Dimension})
        new(@sig(shape), @sig(free_indices), @sig(index_dimensions))
    end
end

is_literal(::Zero) = true
evaluate(::Zero, mapping, component, index_values) = 0

function Base.:(==)(z::Zero, z2::Zero) 
    z.ufl_shape == z2.ufl_shape && z.ufl_free_indices == z2.ufl_free_indices && z.ufl_index_dimensions == z2.ufl_index_dimensions
end
Base.:(==)(::Zero, z2::Real) = z2 == 0 
Base.:(==)(z::Real, z2::Zero) = z2 == z

function Base.repr(z::Zero) 
    if isempty(z.ufl_shape) && isempty(z.ufl_free_indices)
        "0"
    elseif isempty(z.ufl_free_indices)
        "0 shape $(z.ufl_shape)"
    elseif isempty(z.ufl_shape)
        "0 (index labels $(z.ufl_free_indices)"
    else 
        "0 (shape $(z.ufl_shape) index labels $(z.ufl_free_indices)"
    end
end

Base.show(io::IO, ::Zero) = print(io, "0")


@ufl_type struct Identity <: AbstractConstantValue
    ufl_fields = (shape,)
    dim::Dimension 

    function Identity(dim::Dimension)
        new(@sig((dim, dim)), dim)
    end
end

Base.repr(id::Identity) = "Identity($(id.dim))"
Base.show(io::IO, id::Identity) = print(io, "I")
Base.getindex(id::Identity, i::Int, j::Int) = ScalarValue(i == j ? 1 : 0)
Base.getindex(id::Identity, i::Int, j::FixedIndex) = id[i, j.d]
Base.getindex(id::Identity, j::FixedIndex, i::Int) = id[j.d, i]
Base.getindex(id::Identity, i::FixedIndex, j::FixedIndex) = id[i.d, j.d]


@ufl_type struct ScalarValue{T <: Real} <: AbstractConstantValue 
    val::T 

    function ScalarValue(x)
        if x === 0 
            Zero() 
        end

        new{typeof(x)}(@sig(x))
    end
end

Base.show(io::IO, k::ScalarValue) = print(io, k.val)


# Yes this might not handle floating point problems
# But frankly at the moment I could not care 
Base.:(==)(a::ScalarValue{T}, b::T) where T <: Real = a.val == b 
Base.:(==)(a::T, b::ScalarValue{T}) where T <: Real = b == a 


as_ufl(x::AbstractExpr) = x 
as_ufl(x::MultiIndex) = MultiIndexNode(x)
as_ufl(x::Real) = ScalarValue(x)
as_ufl(x) = error("unsupported UFL type for $(x) with type $(typeof(x))")