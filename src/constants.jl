export Zero

abstract type AbstractConstantValue <: Terminal end 

is_cellwise_constant(::AbstractConstantValue) = True 
ufl_domains(::AbstractConstantValue) = ()

struct Zero <: AbstractConstantValue
    @ufl_type 

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