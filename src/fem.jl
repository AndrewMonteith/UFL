export FiniteElement, MixedElement, VectorElement, fem_cell, fem_family, fem_degree, fem_value_shape

abstract type AbstractFiniteElement end 

struct FemBase 
    family::String 
    cell 
    degree::Union{Dimension, DimensionTuple} 
    value_shape::DimensionTuple 
end

fem_family(a::AbstractFiniteElement) = a.base.family 
fem_cell(a::AbstractFiniteElement) = a.base.cell 
fem_degree(a::AbstractFiniteElement) = a.base.degree 
fem_value_shape(a::AbstractFiniteElement) = a.base.value_shape

family_analogies = Dict("CG" => "Lagrange")

family_to_value_ranks = Dict(
    "Lagrange" => 0,
)

struct FiniteElement <: AbstractFiniteElement 
    base::FemBase 

    function FiniteElement(family::String, cell::Cell, degree::Union{Dimension, DimensionTuple})
        family = get(family_analogies, family, family)

        if !(family in keys(family_to_value_ranks))
            error("unknown family $(family)")
        end

        value_rank = family_to_value_ranks[family]
        value_shape = if value_rank === 0
            ()
        elseif value_rank === 1
            (cell.geometric_dimension,)
        elseif value_rank === 2 
            (cell.geometric_dimension, cell.geometric_dimension)
        else
            error("value rank is too high $(value_rank) for $(family). Must be 0 <= r <= 2")
        end

        new(FemBase(family, cell, degree, value_shape))
    end 
end

function maximum_degree(degrees)
    max_len = maximum(length.(degrees))

    degree::Array{Dimension, 1} = [] 
    for ii ∈ 1:max_len 
        push!(degree, 0)
        for j ∈ 1:max_len 
            if length(degrees[j]) >= ii 
                degree[ii] = max(degree[ii], degrees[j][ii])
            end
        end
    end

    if length(degree) === 1
        degree[1]
    else
        tuple(degree...)
    end
end

function NewMixedElementBase(elements...; kwargs...)
    value_shape = get(kwargs, :value_shape, tuple(sum(prod(fem_value_shape(element)) for element ∈ elements)))
    element_family = get(kwargs, :family, "Mixed")

    degree = maximum_degree([fem_degree(element) for element ∈ elements])

    #  In firedrake cell selection is done by sorting them... I'm'a just pick the first one 
    cell = fem_cell(elements[1])

    FemBase(element_family, cell, degree, value_shape)
end

struct MixedElement <: AbstractFiniteElement 
    base::FemBase 

    function MixedElement(elements...; kwargs...)
        new(NewMixedElementBase(elements...; kwargs...))
    end
end

struct VectorElement  <: AbstractFiniteElement
    base::FemBase 

    function VectorElement(element::AbstractFiniteElement; kwargs...)
        cell = fem_cell(element)
        gdim = get(kwargs, :dim, geometric_dimension(cell))

        sub_elements = collect(element for _ ∈ 1:gdim)

        value_shape = tuple(gdim, fem_value_shape(element)...)

        new(NewMixedElementBase(sub_elements...; value_shape=value_shape, family=fem_family(element)))
    end 

    function VectorElement(family::String, cell::Cell, degree::Union{Dimension, DimensionTuple}; kwargs...)
        sub_element = FiniteElement(family, cell, degree)

        VectorElement(sub_element; kwargs...)
    end
end

