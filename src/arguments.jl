export TrialFunction, TestFunction, Constant, UflFunction

abstract type AbstractFormArgument <: Terminal end 

ufl_function_space(arg::AbstractFormArgument) = arg.ufl_function_space
ufl_element(arg::AbstractFormArgument) = arg.ufl_function_space.element
ufl_domain(x::AbstractExpr) = nothing
ufl_domain(arg::AbstractFormArgument) = arg.ufl_function_space.mesh
geometric_dimension(arg::AbstractFormArgument) = ufl_shape(arg)[1]

argument_id = 0
@ufl_type struct UflFunction <: AbstractFormArgument 
    ufl_fields = (shape,)
    id::Int
    ufl_function_space::FunctionSpace

    function UflFunction(fs::FunctionSpace)
        c = argument_id 
        global argument_id = argument_id + 1
        new(fem_value_shape(fs.element), c, fs)
    end
end
Base.show(io::IO, f::UflFunction) = print(io, "w_$(f.id >= 10 ? "{$(f.id)}" : string(f.id))")

@ufl_type struct Argument <: AbstractFormArgument 
    ufl_fields = (shape,)

    number::Int
    part

    ufl_function_space::FunctionSpace
    
    function Argument(function_space::FunctionSpace, number::Int, part = nothing)
        if part !== nothing && !isa(part, Int)
            error("part must be integral or nothing")
        end

        new(@sig(function_space.ufl_shape), @sig(number), @sig(part), function_space)
    end
end
Base.repr(arg::Argument) = "Argument($(repr(arg.ufl_function_space)), $(arg.number))"

is_cellwise_constant(::Argument) = false 

function Base.show(io::IO, arg::Argument)
    s = "v_"
    
    s *= ((arg.number < 10) ? string(arg.number) : "{$(arg.number)}")
    if arg.part !== nothing 
        s *= (arg.part >= 10 ? "^$(arg.part)" : "^{$(arg.part)}")
    end
        
    show(io, s) 
end


function TestFunction(function_space::FunctionSpace, part = nothing)
    Argument(function_space, 0, part)
end 

function TrialFunction(function_space::FunctionSpace, part = nothing)
    Argument(function_space, 1, part)
end

@ufl_type struct Constant <: AbstractFormArgument 
    ufl_fields = (shape,)

    id::Int
    value
    ufl_function_space::FunctionSpace

    function Constant(value::T; @opt(cell::Cell)) where T <: Union{Real, NTuple{N, Real} where N}
        c = argument_id 
        global argument_id = argument_id + 1

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

        new(shape, @sig(c), value, @sig(function_space))
    end
end

Base.repr(c::Constant) = "Constant($((repr ∘ ufl_element)(c.ufl_function_space)), $(c.id)))"
Base.show(io::IO, c::Constant) = print(io, "w_$( c.id >= 10 ? "{$(string(c.id))}" : string(c.id) )")