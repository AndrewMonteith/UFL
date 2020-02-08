ufl_operands(x::Operator)::VarTuple{AbstractExpr} = x.ufl_operands 
ufl_operands(x::Terminal)::VarTuple{AbstractExpr} = ()

struct NoShape end 
struct HasShape end
struct InheritsShape end 

shape_type(::Type{<:Operator}) = NoShape()
shape_type(::Type{ScalarValue{T}}) where T = NoShape()
shape_type(::Type{Indexed}) = NoShape()
shape_type(::Type{MultiIndexNode}) = NoShape()
shape_type(::Type{<:Terminal}) = HasShape()
shape_type(::Type{ComponentTensor}) = HasShape()
shape_type(::Type{IndexSum}) = InheritsShape() 
shape_type(::Type{Sum}) = InheritsShape() 

get_shape(::HasShape, x::AbstractExpr)::DimensionTuple = x.ufl_shape 
get_shape(::NoShape, x::AbstractExpr)::DimensionTuple = ()
get_shape(::InheritsShape, x::AbstractExpr)::DimensionTuple = ufl_shape( ufl_operands(x)[1] )

ufl_shape(x::T) where T = get_shape(shape_type(T), x)
ufl_shape(lt::ListTensor)::DimensionTuple = tuple(length(lt.ufl_operands), ufl_shape(lt.ufl_operands[1])...)


struct HasFreeIndices end 
struct NoFreeIndices end 
struct InheritsFreeIndices end 

free_indices_type(::Type{<:Operator}) = HasFreeIndices()
free_indices_type(::Type{Zero}) = HasFreeIndices()
free_indices_type(::Type{<:Terminal}) = NoFreeIndices() 
free_indices_type(::Type{ListTensor}) = NoFreeIndices()
free_indices_type(::Type{Sum}) = InheritsFreeIndices()
free_indices_type(::Type{Division}) = InheritsFreeIndices()
free_indices_type(::Type{Power}) = InheritsFreeIndices()

get_free_indices(::HasFreeIndices, x::AbstractExpr)::VarTuple{AbstractIndex} = x.ufl_free_indices 
get_free_indices(::NoFreeIndices, x::AbstractExpr)::VarTuple{AbstractIndex} = () 
get_free_indices(::InheritsFreeIndices, x::AbstractExpr) = ufl_free_indices(ufl_operands(x)[1])

ufl_free_indices(x::T) where T = get_free_indices(free_indices_type(T), x)

get_index_dimensions(::HasFreeIndices, x::AbstractExpr)::DimensionTuple = x.ufl_index_dimensions 
get_index_dimensions(::NoFreeIndices, x::AbstractExpr)::DimensionTuple = ()
get_index_dimensions(::InheritsFreeIndices, x::AbstractExpr)::DimensionTuple = ufl_index_dimensions(ufl_operands(x)[1])

ufl_index_dimensions(x::T) where T = get_index_dimensions(free_indices_type(T), x)