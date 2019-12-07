export FixedIndex, Index, MultiIndex

# We do not subtype Operator or Terminal as this will never be directly a node a UFL tree.
abstract type AbstractIndex end 

"""
    An index with a specific value d assigned to it
"""
struct FixedIndex <: AbstractIndex 
    d::Dimension
end 
Base.show(io::IO, i::FixedIndex) = print(io, "FixedIndex($(i.d))")

"""
    Index with a value that ranges between 1 and a dimension d.
    d cannot be initally known hence the data will be provided later 
    alongside the index dimensions
"""
struct Index <: AbstractIndex 
    letter::Char
end 
Base.show(io::IO, i::Index) = print(io, i.letter)

"""
    Sequence of indicies either fixed or free 
    This type is no not a subtype of Terminal 
    This could be a problem? However for know just treating 
    this as a Tuple class seems sufficient.
"""
const MultiIndex = VarTuple{AbstractIndex}

# struct MultiIndex <: Terminal
#     indicies::VarTuple{AbstractIndex}

#     function MultiIndex(indicies::VarTuple{AbstractIndex})
#         new(indicies)
#     end
# end
# Base.length(m::MultiIndex) = Base.length(m.indicies)
# Base.iterate(m::MultiIndex) = Base.iterate(m.indicies)
# Base.show(io::IO, m::MultiIndex) = print(io, "(", join(m.indicies, ", "), ")")
