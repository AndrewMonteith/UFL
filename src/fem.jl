abstract type AbstractFiniteElement end 

struct FemBase 
    family::String 
    cell 
    degree::Union{Dimension, DimensionTuple} 
    value_shape::DimensionTuple 
end

family(a::AbstractFiniteElement) = a.base.family 
cell(a::AbstractFiniteElement) = a.base.cell 
degree(a::AbstractFiniteElement) = a.base.degree 
value_shape(a::AbstractFiniteElement) = a.base.value_shape


struct FiniteElement <: AbstractFiniteElement 
    base::FemBase 

    function FiniteElement(family::String, cell::AbstractCell, degree::Union{Dimension, DimensionTuple})
    end 
end