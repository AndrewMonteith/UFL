export Zero, Identity, IntValue

abstract type AbstractConstantValue <: Terminal end 

is_cellwise_constant(::AbstractConstantValue) = true
ufl_domains(::AbstractConstantValue) = ()

@ufl_type struct Zero <: AbstractConstantValue
    ufl_fields = (shape, free_indices, index_dimensions)

    function Zero(shape::DimensionTuple=())
        new(shape, (), ())
    end

    """
        Constructor for the old format 
        Example old format code:
            i = Index(2)
            j = Index(3)

            z = Zero((), (j, i). Dict{i => 3, j => 5})
    """
    function Zero(shape::DimensionTuple, freeIndices::VarTuple{Index}, indexDimensions::Dict{Index, Dimension})
        sortedIndicies = sort(collect(freeIndices); by=x -> x.id)
        dimensions = Tuple(indexDimensions[index] for index in sortedIndicies)
        new(shape, Tuple(sortedIndicies), dimensions)
    end

    # New Format
    """
        Constructor for the new format 
        Example new format code:
            z = Zero((), (2, 4), (3, 5))
    """
    function Zero(shape::DimensionTuple, freeIndices::VarTuple{Dimension}, indexDimensions::VarTuple{Dimension})
        new(shape, freeIndices, indexDimensions)
    end
end

is_literal(::Zero) = true
evaluate(::Zero, mapping, component, index_values) = 0

function Base.:(==)(z::Zero, z2::Zero) 
    z.ufl_shape == z2.ufl_shape && z.ufl_free_indices == z2.ufl_free_indices && z.ufl_index_dimensions == z2.ufl_index_dimensions
end
Base.:(==)(::Zero, z2::Real) = z2 == 0 
Base.:(==)(z::Real, z2::Zero) = z2 == z

function pretty_print(z::Zero) 
    if z.ufl_shape === () && z.ufl_free_indices === () 
        "0"
    elseif z.ufl_free_indices === ()
        "0 shape $(z.ufl_shape)"
    elseif z.ufl_shape === () 
        "0 (index labels $(z.ufl_free_indices)"
    else 
        "0 (shape $(z.ufl_shape) index labels $(z.ufl_free_indices)"
    end
end

mathsy_print(::Zero) = "0"


@ufl_type struct IntValue <: AbstractConstantValue
    val::Int
end
pretty_print(i::IntValue) = "IntValue($(i.val))"
mathsy_print(i::IntValue) = i.val 

Base.:(==)(i::IntValue, j::Int) = i.val === j 
Base.:(==)(i::Int, j::IntValue) = i === j.val 

@ufl_type struct Identity <: AbstractConstantValue
    ufl_fields = (shape,)

    dim::Dimension 

    function Identity(dim::Dimension)
        new(dim, (dim, dim))
    end
end

Base.show(i::IO, id::Identity) = print("Identity($(id.dim))")
Base.getindex(id::Identity, i::Int, j::Int) = IntValue(i == j ? 1 : 0)
Base.getindex(id::Identity, i::Int, j::FixedIndex) = id[i, j.d]
Base.getindex(id::Identity, j::FixedIndex, i::Int) = id[j.d, i]
Base.getindex(id::Identity, i::FixedIndex, j::FixedIndex) = id[i.d, j.d]