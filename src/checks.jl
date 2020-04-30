export is_true_ufl_scalar # , is_scalar_constant_expression

# is_scalar_constant_expression(r::Real) = true
# is_scalar_constant_expression(e::AbstractConstantValue) = true 
# is_scalar_constant_expression(e::AbstractFormArgument) = ufl_element(e).family === "Real"
# function is_scalar_constant_expression(e::AbstractExpr)
#     !isempty(e) && return false 

#     UFL.@pre_order_traversal for op ∈ e 
#         if !(is_scalar_constant_expression)
#             (op <: Operator || op <: AbstractConstantValue) && continue 

#             if op <: AbstractFormArgument 
#                 ufl_element(e).family === "Real" || return false 
#             end
#         end

#         return true 
#     end

is_true_ufl_scalar(e::AbstractExpr) = (isempty ∘ ufl_shape)(e) && (isempty ∘ ufl_free_indices)(e)