export FiniteElement, MixedElement, VectorElement, fem_cell, fem_family, fem_degree, fem_value_shape, fem_mapping, fem_ref_value_shape

abstract type AbstractFiniteElement end 

struct FemBase 
    family::String 
    cell::UFL.@opt_t(Cell)
    degree::UFL.@opt_t(Union{Dimension, DimensionTuple})
    value_shape::DimensionTuple 
    ref_value_shape::DimensionTuple
    mapping::String

    function FemBase(family::String, cell::UFL.@opt_t(Cell), degree::UFL.@opt_t(Union{Dimension, DimensionTuple}), value_shape::DimensionTuple, ref_val_shape::DimensionTuple)
        mapping = get!(family_to_mappings, family, "affine")
        new(family, cell, degree, value_shape, ref_val_shape, mapping)
    end
end

family_analogies = Dict("CG" => "Lagrange", "R" => "Real")

family_to_value_ranks = Dict(
    "Lagrange" => 0,
    "Real" => 0
)

family_to_mappings = Dict(
    "Lagrange" => "identity",
    "Real" => "identity"
)


struct FiniteElement <: AbstractFiniteElement 
    base::FemBase 

    function FiniteElement(family::String; @opt(cell::Cell), @opt(degree::Union{Dimension, DimensionTuple}))
        family = get(family_analogies, family, family)

        if !(family in keys(family_to_value_ranks))
            error("unknown family $(family)")
        end

        value_rank = family_to_value_ranks[family]
        value_shape, ref_value_shape = if value_rank === 0
            (), ()
        elseif cell === nothing 
            error("cannot infer shape with no provided cell")
        elseif value_rank === 1
            (cell.geometric_dimension,), (cell.topological_dimension,)
        elseif value_rank === 2 
            (cell.geometric_dimension, cell.geometric_dimension), (cell.topological_dimension, cell.topological_dimension)
        else
            error("value rank is too high $(value_rank) for $(family). Must be 0 <= r <= 2")
        end

        new(FemBase(family, cell, degree, value_shape, ref_value_shape))
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
    value_shape = get(kwargs, :value_shape, tuple(sum(collect(prod(fem_value_shape(e)) for e ∈ elements)))) 
    ref_value_shape = get(kwargs, :ref_value_shape, tuple(sum(collect(prod(fem_ref_value_shape(e)) for e ∈ elements)))) 
    element_family = get(kwargs, :family, "Mixed")

    degree = maximum_degree([fem_degree(element) for element ∈ elements])

    #  In firedrake cell selection is done by sorting them... I'm'a just pick the first one 
    cell = fem_cell(elements[1])

    FemBase(element_family, cell, degree, value_shape, ref_value_shape)
end

struct MixedElement <: AbstractFiniteElement 
    base::FemBase 
    elements::VarTuple{AbstractFiniteElement}

    function MixedElement(elements...; kwargs...)
        new(NewMixedElementBase(elements...; kwargs...), elements)
    end
end

struct VectorElement  <: AbstractFiniteElement
    base::FemBase 

    function VectorElement(family::String; @opt(cell::Cell), @opt(degree::Union{Dimension, DimensionTuple}), @opt(dim::Dimension))
        sub_element = FiniteElement(family; cell=cell, degree=degree)

        VectorElement(sub_element; dim=dim)
    end

    function VectorElement(element::AbstractFiniteElement; @opt(dim::Dimension))
        cell = fem_cell(element)

        if dim === nothing 
            dim = geometric_dimension(cell)
        end

        sub_elements = collect(element for _ ∈ 1:dim)

        value_shape = tuple(dim, fem_value_shape(element)...)

        new(NewMixedElementBase(sub_elements...; value_shape=value_shape, family=fem_family(element)))
    end 
end


fem_family(a::AbstractFiniteElement) = a.base.family 
fem_cell(a::AbstractFiniteElement) = a.base.cell 
fem_degree(a::AbstractFiniteElement) = a.base.degree 
fem_value_shape(a::AbstractFiniteElement) = a.base.value_shape
fem_ref_value_shape(a::AbstractFiniteElement) = a.base.ref_value_shape
fem_mapping(a::AbstractFiniteElement) = a.base.mapping 
fem_mapping(m::MixedElement) = all(fem_mapping(e)==="identity" for e ∈ m.elements) ? "identity" : "undefined"