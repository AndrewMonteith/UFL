export AbstractIndex, FixedIndex, Index, MultiIndex

# We do not subtype Operator or Terminal as this will never be directly a node a UFL tree.
abstract type AbstractIndex end 

"""
    An index with a specific value d assigned to it
"""
struct FixedIndex <: AbstractIndex 
    d::Dimension
end 
Base.show(io::IO, i::FixedIndex) = print(io, "FixedIndex($(i.d))")
# Base.convert(::Type{AbstractIndex}, i::FixedIndex) = i

"""
    An index with no value assigned ie can represent some value between [1, d]
    where d is a dimension
"""
struct Index <: AbstractIndex 
    d::Dimension
end 
Base.show(io::IO, i::Index) = print(io, "Index($(i.d))")

"""
    Sequence of indicies either fixed or free 

    N.B: Thought about having the member variable as a Tuple.
         Julias tuples are not homogenous in the type of each element so it's included in the type parameters
         This makes it unsuitable as we don't know the number of dimensions, Julia provides the NTuple 
         class which are homogeonous but you need to provide it a size which means I would need 
         parameterize MultiIndex. For this reason I chose to use an array, whilst they're mutable 
         I don't foresee this being a problem in the future.
"""
struct MultiIndex <: Terminal 
    indicies::Array{AbstractIndex}

    function MultiIndex(indicies::NTuple{N, AbstractIndex}) where N 
        new(collect(indicies))
    end

end
Base.show(io::IO, m::MultiIndex) = print(io, "(", join(m.indicies, ", "), ")")
Base.length(m::MultiIndex) = length(m.indicies)
Base.iterate(m::MultiIndex) = iterate(m.indicies)
# convert(::Type{MultiIndex}, indicies::NTuple{N, AbstractIndex}) where N = MultiIndex(indicies)
convert(::Type{MultiIndex}, indicies::Tuple{Index}) = MultiIndex(indicies)

ufl_shape(::MultiIndex) = error("MultiIndex has no shape")
ufl_free_indicies(::MultiIndex) = error("MultiIndex has no free indicies")
ufl_index_dimensions(::MultiIndex) = error("MultiIndex has no free indicies")
