export ufl_operands 

abstract type Operands end
struct HasOperands <: Operands end 
struct NoOperands <: Operands end 

has_operands(::Type{<:Operator}) = HasOperands()
has_operands(::Type{<:Terminal}) = NoOperands()

ufl_operands(x::T) where T = ufl_operands(has_operands(T), x)
ufl_operands(::HasOperands, x) = x.ufl_operands 
ufl_operands(::NoOperands, x) = ()

# Is it worth going down this route more?