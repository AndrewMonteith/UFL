export Zero

abstract type AbstractConstantValue <: Terminal end 

is_cellwise_constant(::AbstractConstantValue) = True 
ufl_domains(::AbstractConstantValue) = ()

struct Zero <: AbstractConstantValue
    @ufl_type(Zero, ufl_shape, ufl_free_indices, ufl_index_dimensions) 

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

function Base.show(io::IO, z::Zero) 
    if z.ufl_shape === () && self.ufl_free_indices === () 
        print(io, "0")
    elseif z.ufl_free_indices === ()
        print(io, "0 shape $(z.ufl_shape)")
    elseif z.ufl_shape === () 
        print(io, "0 (index labels $(z.ufl_free_indices)")
    else 
        print(io, "0 (shape $(z.ufl_shape) index labels $(z.ufl_free_indices)")
    end
end
