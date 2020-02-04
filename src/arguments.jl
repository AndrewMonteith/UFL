export TrialFunction, TestFunction, Constant

abstract type AbstractFormArgument end 

ufl_function_space(arg::AbstractFormArgument) = arg.ufl_function_space
ufl_domain(arg::AbstractFormArgument) = ufl_element(arg.ufl_function_space)
ufl_domains(arg::AbstractFormArgument) = ufl_domains(arg.ufl_function_space)
geometric_dimension(arg::AbstractFormArgument) = ufl_shape(arg)[1]

@ufl_type struct Argument <: AbstractFormArgument 
    number::Int
    part

    ufl_function_space::FunctionSpace
    
    function Argument(function_space::FunctionSpace, number::Int, part = nothing)
        if part !== nothing && !isa(part, Int)
            error("part must be integral or nothing")
        end

        new(function_space.ufl_shape, (), (), number, part, function_space)
    end
end
Base.repr(arg::Argument) = "Argument($(repr(arg.ufl_function_space)), $(arg.number))"

is_cellwise_constant(::Argument) = false 

function Base.show(io::IO, arg::Argument)
    s = "v_"
    
    s *= arg.number < 10 ? arg.number : "{$(arg.number)}"    
    s *= arg.part !== nothing ? "^$(arg.part)" : "^{$(arg.part)}"
    
    show(io, s) 
end


function TestFunction(function_space::FunctionSpace, part = nothing)
    Argument(function_space, 0, part)
end 

function TrialFunction(function_space::FunctionSpace, part = nothing)
    Argument(function_space, 1, part)
end



@ufl_type struct Constant <: AbstractFormArgument 
    value
    ufl_function_space::FunctionSpace

    function Constant(value::T; @opt(cell::Cell)) where T <: Union{Real, NTuple{N, Real} where N}
        shape = isa(value, Tuple) ? length(value) : size(value)
        rank = length(shape)

        e = if rank === 0
            FiniteElement("Real"; cell=cell, degree=0)
        elseif rank === 1
            VectorElement("Real"; cell=cell, degree=0, dim=shape[1])
        else 
            error("not supporting tensors yet")
        end

        mesh = cell !== nothing ? Mesh(cell) : cell
        function_space = FunctionSpace(e; mesh=mesh)

        shape = fem_value_shape(function_space.element) 

        new(shape, (), (), value, function_space)
    end
end

Base.repr(c::Constant) = "Constant($((repr ∘ ufl_element)(c))))"
Base.show(io::IO, c::Constant) = print(io, repr(c))