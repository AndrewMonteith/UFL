export grad

abstract type AbstractDifferential <: Operator end 

is_differential(::AbstractDifferential) = true 

abstract type CompoundDerivative <: AbstractDifferential end 

#=
    Incomplete type 
    TODO: 
        1. Implement dimension stuff with Grad
            - Current Block: No primitives that have any dimension 
=#
@ufl_type struct Grad <: CompoundDerivative 
    ufl_fields = (operands,)
    ufl_tags = (inherit_indices_from_operand=1,)

    function Grad(expr::AbstractExpr, operands::VarTuple{AbstractExpr})
        # TODO: Simplification if expr is_cellwise_constant
        #   - Requirements: find_geometric_dimension (ie dimension stuff)
        #                   Could test whether just directly calling geometric_dimension works

        new(operands)
    end
end

num_of_ops(::Grad) = 1
ufl_shape(g::Grad) = tuple(ufl_shape(ufl_operands(g)[1])..., g.dim)

function grad(f)
    Grad(as_ufl(f))
end