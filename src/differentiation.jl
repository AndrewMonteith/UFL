export grad, derivative, Coefficient

abstract type AbstractDerivative <: Operator end 

const Coefficient = Union{UflFunction,Constant}

@ufl_type struct Grad <: AbstractDerivative
    ufl_fields = (operands, shape)
    ufl_tags = (num_ops = 1,)

    dim::Dimension

    function Grad(f::AbstractExpr)
        is_cellwise_constant(f) && return Zero(tuple(ufl_shape(f)..., geometric_dimension(f)), ufl_free_indices(f), ufl_index_dimensions(f))

        gdim = find_geometric_dimension(f)
        new(@sig((f,)), tuple(ufl_shape(f)..., gdim), gdim)
    end
end
Base.show(io::IO, g::Grad) = print(io, "grad($(g.ufl_operands[1]))")

grad(f)::AbstractExpr = (Grad ∘ as_ufl)(f)
function reconstruct_expr(g::Grad, op::AbstractExpr)::AbstractExpr
    if is_cellwise_constant(op)
        ufl_shape(op) !== ufl_shape(g.ufl_operands[1]) && error("Operands shape mismatch in Grad reconstruct")
        ufl_free_indices(op) !== ufl_free_indices(g.ufl_operands[1]) && error("Free index mismatch in Grad reconstruct")

        Zero(ufl_shape(g), ufl_free_indices(g), ufl_index_dimensions(g))
    else
        grad(op)
    end 
end 


@ufl_type struct ReferenceGrad <: AbstractDerivative
    ufl_fields = (operands,)
    dim::Dimension 

    function ReferenceGrad(f::AbstractExpr)
        dim = (topological_dimension ∘ ufl_domain)(f)

        if is_cellwise_constant(f)
            return Zero(tuple(f.ufl_shape..., dim), ufl_free_indices(f), ufl_index_dimensions(f))
        end

        new(@sig((f,)), dim)
    end
end 
ufl_shape(rg::ReferenceGrad) = tuple(ufl_shape(rg.ufl_operands[1])..., rg.dim)


@ufl_type struct CoefficientDerivative <: AbstractDerivative
    ufl_fields = (operands,)
    ufl_tags = (num_ops = 4,)

    function CoefficientDerivative(integrand::AbstractExpr, coefficients::ExprList, arguments::ExprList, coefficient_derivatves::ExprList)
        new(@sig((integrand, coefficients, arguments, coefficient_derivatves)))
    end 
end 
function Base.show(io::IO, cd::CoefficientDerivative)
    a, b, c, d = cd.ufl_operands 
    s = "d/dfj { $a }, with fh=$b, dfh/dfj = $c and coefficient derivatives $d"

    print(io, s)
end

function handle_derivative_arguments(form::Union{Form,AbstractExpr}, coefficient::Coefficient; @opt(argument::Argument))
    # This function has been heavily simplified and may need further development in the future.
    # Simplied for the case we have 1 coefficient and 1 argument 

    ufl_shape(coefficient) !== ufl_shape(argument) && error("coefficient and argument cannot have mismatching shapes")

    coefficients = ExprList(coefficient)
    arguments = ExprList(argument)

    return coefficients, arguments
end

function gateaux_derivative(form::Union{Form,AbstractExpr}, coefficient::Coefficient; @opt(argument::Argument), @opt(coefficient_derivatives::Dict{Coefficient,Argument}))
    # Compute Gateaux derivative for form w.r.t coefficient in direction of argument 
    coefficients, arguments = handle_derivative_arguments(form, coefficient; argument = argument)

    if coefficient_derivatives === nothing 
        coefficient_derivatives = ExprList()
    else
        error("need to implement this")
    end 

    if form isa Form 
        integrals = [] 
        for integral ∈ form.integrals 
            # TODO? : Accept coefficient as SpatialCoordinate 
            fd = CoefficientDerivative(integral.integrand, coefficients, arguments, coefficient_derivatives)
            push!(integrals, reconstruct(integral; integrand = fd))
        end 

        Form(tuple(integrals...))
    elseif form isa AbstractExpr 
        CoefficientDerivative(form, coefficients, arguments, coefficient_derivatives)
    else
        error("Invalid argument type")
    end 
end

function derivative(form::Form, u::Coefficient; @opt(du::Argument))
    function argument(V::FunctionSpace)
        if du === nothing 
            n = isempty(form.arguments) ? -1 : maximum(u.number for arg ∈ form.arguments) 
            Argument(V, n + 1)
        else
            du 
        end 
    end

    du = if u isa UflFunction 
        (argument ∘ ufl_function_space)(u)
    elseif u isa Constant 
        !(isempty ∘ ufl_shape)(u) && error("Real function space of vector elements not supported")

        V = FunctionSpace(form.domain, "Real", 0) 
        argument(V)
    end 
    
    gateaux_derivative(form, u; argument = du)
end 