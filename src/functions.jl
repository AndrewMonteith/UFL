# geometric_dimension would be exported by dimensions.jl if we had one
export Argument, TrialFunction, TestFunction, FunctionSpace, geometric_dimension

"""
    Mishmash of datatypes needed to get TrialFunction & TestFunction 
    working and bits and bobs. Eventually this file will get separated
    into the appropriate files when I work out what FEM stuff 
"""
abstract type AbstractFormArgument <: Terminal end 


@ufl_type struct Argument <: AbstractFormArgument 
    ufl_fields = (shape, function_space)
    
    number::Int
    part
    
    function Argument(function_space, number::Int, part = nothing)
        if part !== nothing && !isa(part, Int)
            error("part must be integral or nothing")
        end


        # function_space should be a subtype of AbstractFunctionSpace
        new(number, part, ufl_shape(function_space), function_space)
    end
end

ufl_domain(arg::Argument) = ufl_element(arg.ufl_function_space)
ufl_domains(arg::Argument) = ufl_domains(arg.ufl_function_space)
geometric_dimension(arg::Argument) = ufl_shape(arg)[1]
is_cellwise_constant(::Argument) = false 

Base.repr(io::IO, arg::Argument) = repr(io, "Argument(shape $(ufl_shape(arg)) number $(arg.number) part $(arg.part)")
function Base.show(io::IO, arg::Argument)
    s = "v_"
    
    s *= arg.number < 10 ? arg.number : "{$(arg.number)}"    
    s *= arg.part !== nothing ? "^$(arg.part)" : "^{$(arg.part)}"
    
    show(io, s) 
end


abstract type AbstractFunctionSpace end 

# Dummy function space until I can work out what I actually need 
struct FunctionSpace <: AbstractFunctionSpace end

ufl_domain(::FunctionSpace) = error("not implemented yet")
ufl_domains(::FunctionSpace) = error("not implemented yet")
ufl_shape(::FunctionSpace) = (2,)


function TestFunction(function_space::AbstractFunctionSpace, part = nothing)
    Argument(function_space, 0, part)
end 

function TrialFunction(function_space::AbstractFunctionSpace, part = nothing)
    Argument(function_space, 1, part)
end


