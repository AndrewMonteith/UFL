export Terminal, Operator, Dimension, geometric_dimension, topological_dimension, compute_hash, ufl_shape, ufl_operands, ufl_free_indices, ufl_index_dimensions, hash_behaviour

"""
    Root type of any node in the UFL tree.
"""
abstract type AbstractExpr end 

"""
An expression that does not depend on any other expression.
Typically a Terminal has some non-expression data associated with it 
such as geometry data or constants.
"""
abstract type Terminal <: AbstractExpr end 

"""
    A result of an operator, such as IndexSum, ComponentTensor, MinValue, ...
    If Terminal types represnet leaf nodes of the type tree, then Operator types
    will represent the non-terminal types.
"""
abstract type Operator <: AbstractExpr end 

# compute_hash(o::Operator) = hash((ufl_typecode(o), (hash(op) for op in ufl_operands(o))...))
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


function ufl_shape end 
function ufl_operands end
function ufl_free_indices end 
function ufl_index_dimensions end 
function reconstruct_expr end

function reconstruct_expr(t::Terminal, operands::VarTuple{AbstractExpr})
    !isempty(operands) && error("Terminal has no operands")
    t
end 

function reconstruct_expr(o::Operator, operands::VarTuple{AbstractExpr})
    typeof(o)(operands...)
end

is_cellwise_constant(::AbstractExpr) = false
geometric_dimension(::AbstractExpr) = -1

function parstr end 

hash_behaviour(x::Any) = x 
hash_behaviour(x::AbstractExpr) = x.ufl_hash_code

topological_dimension(x::Any)::Dimension = x.topological_dimension