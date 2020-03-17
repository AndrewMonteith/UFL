# struct BinaryOperator end 
# struct HasOperands end 
# struct NoOperands end 

# operands_type(::Type{<:Operator}) = HasOperands()
# operands_type(::Type{<:Terminal}) = NoOperands() 
# operands_type(::Type{Sum}) = BinaryOperator()
# operands_type(::Type{Product}) = BinaryOperator()

# Type Instability of these functions makes the program a bit slower:
# get_operands(::HasOperands, x::AbstractExpr)::VarTuple{AbstractExpr} = x.ufl_operands 
# get_operands(::NoOperands, x::AbstractExpr)::VarTuple{AbstractExpr} = () 
# get_operands(::BinaryOperator, x::AbstractExpr)::Tuple{AbstractExpr, AbstractExpr} = x.ufl_operands

# ufl_operands(x::T) where T <: AbstractExpr = get_operands(operands_type(T), x)

ufl_operands(x::Operator)::VarTuple{AbstractExpr} = x.ufl_operands 
ufl_operands(x::Terminal)::VarTuple{AbstractExpr} = ()

struct NoShape end 
struct HasShape end
struct InheritsShape end 

shape_type(::Type{<:Operator}) = NoShape()
shape_type(::Type{ScalarValue{T}}) where T = NoShape()
shape_type(::Type{Indexed}) = NoShape()
shape_type(::Type{MultiIndexNode}) = NoShape()
shape_type(::Type{Trace}) = NoShape()
shape_type(::Type{Determinant}) = NoShape()
shape_type(::Type{<:Terminal}) = HasShape()
shape_type(::Type{ComponentTensor}) = HasShape()
shape_type(::Type{Grad}) = HasShape()
shape_type(::Type{UFL.ReferenceValue}) = HasShape()
shape_type(::Type{<:CompoundTensorOperator}) = HasShape()
shape_type(::Type{Jacobian}) = HasShape()

shape_type(::Type{IndexSum}) = InheritsShape() 
shape_type(::Type{Sum}) = InheritsShape() 
shape_type(::Type{CoefficientDerivative}) = InheritsShape()
shape_type(::Type{ExprList}) = error("non-tensor type has no shape")

get_shape(::HasShape, x::AbstractExpr)::DimensionTuple = x.ufl_shape 
get_shape(::NoShape, x::AbstractExpr)::DimensionTuple = ()
get_shape(::InheritsShape, x::AbstractExpr)::DimensionTuple = ufl_shape( ufl_operands(x)[1] )

ufl_shape(x::T) where T = get_shape(shape_type(T), x)
ufl_shape(lt::ListTensor)::DimensionTuple = tuple(length(lt.ufl_operands), ufl_shape(lt.ufl_operands[1])...)


struct HasFreeIndices end 
struct NoFreeIndices end 
struct InheritsFreeIndices end 

free_indices_type(::Type{<:Terminal}) = NoFreeIndices() 
free_indices_type(::Type{ListTensor}) = NoFreeIndices()
free_indices_type(::Type{<:CompoundTensorOperator}) = NoFreeIndices()
free_indices_type(::Type{UFL.ReferenceValue}) = NoFreeIndices()
free_indices_type(::Type{<:MathFunction}) = NoFreeIndices()
free_indices_type(::Type{<:Operator}) = HasFreeIndices()
free_indices_type(::Type{Zero}) = HasFreeIndices()
free_indices_type(::Type{Sum}) = InheritsFreeIndices()
free_indices_type(::Type{Division}) = InheritsFreeIndices()
free_indices_type(::Type{Power}) = InheritsFreeIndices()
free_indices_type(::Type{Grad}) = InheritsFreeIndices()
free_indices_type(::Type{Transposed}) = InheritsFreeIndices()
free_indices_type(::Type{Trace}) = InheritsFreeIndices()
free_indices_type(::Type{CoefficientDerivative}) = InheritsFreeIndices()
free_indices_type(::Type{ExprList}) = error("non-tensor type has no free indicies")

get_free_indices(::HasFreeIndices, x::AbstractExpr)::VarTuple{AbstractIndex} = x.ufl_free_indices 
get_free_indices(::NoFreeIndices, x::AbstractExpr)::VarTuple{AbstractIndex} = () 
get_free_indices(::InheritsFreeIndices, x::AbstractExpr) = ufl_free_indices(ufl_operands(x)[1])

ufl_free_indices(x::T) where T = get_free_indices(free_indices_type(T), x)

get_index_dimensions(::HasFreeIndices, x::AbstractExpr)::DimensionTuple = x.ufl_index_dimensions 
get_index_dimensions(::NoFreeIndices, x::AbstractExpr)::DimensionTuple = ()
get_index_dimensions(::InheritsFreeIndices, x::AbstractExpr)::DimensionTuple = ufl_index_dimensions(ufl_operands(x)[1])

ufl_index_dimensions(x::T) where T = get_index_dimensions(free_indices_type(T), x)


Base.hash(x::AbstractExpr) = x.ufl_hash_code
Base.:(==)(x1::AbstractExpr, x2::AbstractExpr) = hash(x1) === hash(x2)

function Base.length(x::AbstractExpr)
    sh = ufl_shape(x)
    length(sh) > 1 && error("cannot take length of non-vector expression")
    sh[1]
end

function Base.iterate(x::AbstractExpr, state::Tuple{Int, Int})
    state[1] > state[2] && (nothing, nothing)

    i = state[1]+1

    (x[i], (i, state[2]))
end
Base.iterate(x::AbstractExpr) = (x[1], (1, length(x)))