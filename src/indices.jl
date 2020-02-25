export FixedIndex, Index, MultiIndex

# We do not subtype Operator or Terminal as this will never be directly a node a UFL tree.
abstract type AbstractIndex end 

# We apply the ufl_type macros to the indices not because they are part the UFL type tree
# but because we want to apply the hashing behaviour

# An index with a specific value d assigned to it
@ufl_type struct FixedIndex <: AbstractIndex 
    d::Dimension

    function FixedIndex(d::Dimension)
        new(@sig(d))
    end
end
Base.show(io::IO, i::FixedIndex) = print(io, i.d)
Base.repr(i::FixedIndex) = "FixedIndex($(i.d))"

# Index with a value that ranges between 1 and a dimension d.
# d cannot be initally known hence the data will be provided later 
#         alongside the index dimensions
index_count = 0
@ufl_type struct Index <: AbstractIndex 
    id::Int

    function Index()
        id = index_count
        global index_count = index_count + 1

        new(@sig(id))
    end

    """
        Used for testing purposes
    """
    function Index(id::Int)
        new(@sig(id))
    end
end 
Base.show(io::IO, i::Index) = print(io, "i_$(i.id)")
Base.:(==)(i::Index, j::Index) = i.id === j.id
Base.:(==)(i::Index, j::Int) = i.id === j
Base.:(==)(i::Int, j::Index) = j == i
Base.isless(i::Index, j::Index) = i.id < j.id

"""
    Note:
    Sequence of indices either fixed or free 
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
function Base.show(io::IO, m::MultiIndexNode)
    if length(m.indices) === 1
        print(io, m.indices[1])
    else
        print(io, join(m.indices, ", "))
    end 
end
Base.repr(m::MultiIndexNode) = "MultiIndex$(repr(m.indices))"

Base.iterate(m::MultiIndexNode) = isempty(m.indices) ? nothing : (m.indices[1], 1)
Base.iterate(m::MultiIndexNode, state::Int) = state === length(m.indices) ? nothing : (m[state+1], state+1)
convert(::Type{VarTuple{AbstractIndex}}, x) = MultiIndexNode(x)

indices_n(n::Int) = tuple((Index() for _ ∈ 1:n)...)
indices(m::MultiIndexNode) = m.indices
indices(m) = m
