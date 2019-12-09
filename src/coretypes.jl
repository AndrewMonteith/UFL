export Terminal, Operator, Dimension

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