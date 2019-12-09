abstract type AbstractDifferential <: Operator end 

is_differential(::AbstractDifferential) = true 

abstract type CompoundDerivative <: AbstractDifferential end 

@ufl_type struct Grad <: CompoundDerivative 
    ufl_fields = (operands,)

    ufl_tags = (inherit_indices_from_operand=0,)
end

num_of_ops(::Grad) = 1
