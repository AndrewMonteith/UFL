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
    Note to future self:
        May need to maintain a cache of the indices 
        Some places we refer only to the id 
        We may need to recover the original struct if ever
        more information was associated with it 
"""


"""
    Index with a value that ranges between 1 and a dimension d.
    d cannot be initally known hence the data will be provided later 
    alongside the index dimensions
"""
index_count = 0
struct Index <: AbstractIndex 
    id::Int

    function Index()
        id = index_count++

        new(id)
    end

    """
        Used for testing purposes
    """
    function Index(id::Int)
        new(id)
    end
end 
Base.show(io::IO, i::Index) = print(io, "i_$(i.id)")
Base.:(==)(i::Index, j::Index) = i.id === j.id
Base.:(==)(i::Index, j::Int) = i.id === j
Base.:(==)(i::Int, j::Index) = j == i

"""
    Sequence of indices either fixed or free 
    Note:
        This type is not a subtype of Terminal 
        This could be a problem? However for know just treating 
        this as a class seems sufficient.
"""
const MultiIndex = VarTuple{AbstractIndex}

# struct MultiIndex <: Terminal
#     indices::VarTuple{AbstractIndex}

#     function MultiIndex(indices::VarTuple{AbstractIndex})
#         new(indices)
#     end
# end
# Base.length(m::MultiIndex) = Base.length(m.indices)
# Base.iterate(m::MultiIndex) = Base.iterate(m.indices)
# Base.show(io::IO, m::MultiIndex) = print(io, "(", join(m.indices, ", "), ")")
