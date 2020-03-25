export apply_function_pullback

@ufl_type struct ReferenceValue <: Operator
    ufl_fields = (operands, shape)

    function ReferenceValue(f::AbstractFormArgument)
        new(@sig((f,)), (fem_ref_value_shape ∘ ufl_element)(f))
    end
end
Base.show(io::IO, rv::ReferenceValue) = print(io, "reference_value($(rv.ufl_operands[1]))")

function apply_single_function_pullback(arg::AbstractFormArgument)
    element = ufl_element(arg)
    mapping = fem_mapping(element)

    r = ReferenceValue(arg)

    argsh, rsh = ufl_shape(arg), ufl_shape(r)

    if mapping === "identity"
        argsh != rsh && error("this should not happen")
        return r 
    end

    error("undefined for $arg with mapping $mapping")
end 

struct FunctionPullback <: AbstractMapper 
    base::BaseMapper 
    cache::Dict{AbstractFormArgument, AbstractExpr} 

    FunctionPullback() = new(BaseMapper(), Dict{AbstractFormArgument, AbstractExpr}())
end 

(fs::FunctionPullback)(t::Terminal) = t
(fs::FunctionPullback)(expr::AbstractExpr) = fs.base(expr)
(fs::FunctionPullback)(arg::AbstractFormArgument) = get!(fs.cache, arg, apply_single_function_pullback(arg))

apply_function_pullback(expr::Union{Form, AbstractExpr}) = map_integrand_dags(FunctionPullback(), expr)