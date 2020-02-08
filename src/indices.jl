export FixedIndex, Index, MultiIndex

# We do not subtype Operator or Terminal as this will never be directly a node a UFL tree.
abstract type AbstractIndex end 

"""
    An index with a specific value d assigned to it
"""
struct FixedIndex <: AbstractIndex 
    d::Dimension
end 
Base.show(io::IO, i::FixedIndex) = show(io, "FixedIndex($(i.d))")


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
        id = index_count
        global index_count = index_count + 1

        new(id)
    end

    """
        Used for testing purposes
    """
    function Index(id::Int)
        new(id)
    end
end 
Base.show(io::IO, i::Index) = show(io, "i_$(i.id)")
Base.:(==)(i::Index, j::Index) = i.id === j.id
Base.:(==)(i::Index, j::Int) = i.id === j
Base.:(==)(i::Int, j::Index) = j == i
Base.isless(i::Index, j::Index) = i.id < j.id

"""
    Sequence of indices either fixed or free 
    Note:
        This type is not a subtype of Terminal 
        This could be a problem? However for know just treating 
        this as a class seems sufficient. We can just treat this as an 
        edge case
"""
const MultiIndex = VarTuple{AbstractIndex}


@ufl_type struct MultiIndexNode <: Terminal
    indices::VarTuple{AbstractIndex}
    
    function MultiIndexNode(indices::VarTuple{AbstractIndex})
        new(indices)
    end
end

Base.length(m::MultiIndexNode) = Base.length(m.indices)
Base.iterate(m::MultiIndexNode) = Base.iterate(m.indices)
Base.show(io::IO, m::MultiIndexNode) = print(io, "(", join(m.indices, ", "), ")")
Base.repr(m::MultiIndexNode) = "MultiIndex$(repr(m.indices))"

convert(::Type{VarTuple{AbstractIndex}}, x) = MultiIndexNode(x)

indices_n(n::Int) = tuple((Index() for _ ∈ 1:n)...)
indices(m::MultiIndexNode) = m.indices
indices(m) = m
