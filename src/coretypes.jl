export Terminal, Operator, Dimension, geometric_dimension, topological_dimension, ufl_compute_hash

"""
    Root type of any node in the UFL tree.
"""
abstract type AbstractExpr end 

Base.hash(a::AbstractExpr) = compute_expr_hash(a)

"""
An expression that does not depend on any other expression.
Typically a Terminal has some non-expression data associated with it 
such as geometry data or constants.
"""
abstract type Terminal <: AbstractExpr end 

ufl_compute_hash(t::Terminal) = (hash ∘ repr)(t)

"""
    A result of an operator, such as IndexSum, ComponentTensor, MinValue, ...
    If Terminal types represnet leaf nodes of the type tree, then Operator types
    will represent the non-terminal types.
"""
abstract type Operator <: AbstractExpr end 

ufl_compute_hash(o::Operator) = hash((ufl_typecode(o), (hash(op) for op in ufl_operands(o))...))
function Base.repr(o::Operator)
    str_ops = join((repr(op) for op in ufl_operands(o)), ",")
    "$(typeof(o))$(str_ops)"
end

"""
    Variable length tuple of a type T 
"""
const VarTuple{T} = NTuple{N,T} where N


"""
    A dimension represents any strictly positive integer
    We make it Int32 to just shut up any conversion errors
"""
const Dimension = Int32


const DimensionTuple = VarTuple{Dimension}



geometric_dimension(x::Any)::Dimension = x.geometric_dimension
topological_dimension(x::Any)::Dimension = x.topological_dimension