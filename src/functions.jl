# geometric_dimension would be exported by dimensions.jl if we had one
export Argument, TrialFunction, TestFunction, FunctionSpace, geometric_dimension, VectorFunctionSpace



#=
Add required data as needed.
=#
struct FunctionSpace 
    ufl_shape::DimensionTuple
    mesh::Mesh 
    element::AbstractFiniteElement
    
    # degree::Dimension Not sure how to set this... 
    
    function FunctionSpace(mesh::Mesh, element::AbstractFiniteElement)
        shape = if isa(element, MixedElement) 
            error("cannot create FunctionSpace from MixedElement")
        elseif isa(element, VectorElement) 
            fem_value_shape(element)
        else
            tuple(fem_value_shape(element)[:1])
        end 
        
        new(shape, mesh, element)
    end
end

function VectorFunctionSpace(mesh::Mesh, element::AbstractFiniteElement; kwargs...)
    dim = get(kwargs, :dim, geometric_dimension(mesh))

    element = VectorElement(element; dim=dim)

    FunctionSpace(mesh, element)
end


@ufl_type struct Argument <: Terminal 
    ufl_fields = (shape,)
    
    number::Int
    part

    ufl_function_space::FunctionSpace
    
    function Argument(function_space::FunctionSpace, number::Int, part = nothing)
        if part !== nothing && !isa(part, Int)
            error("part must be integral or nothing")
        end

        new(function_space.ufl_shape, number, part, function_space)
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

function TestFunction(function_space::FunctionSpace, part = nothing)
    Argument(function_space, 0, part)
end 

function TrialFunction(function_space::FunctionSpace, part = nothing)
    Argument(function_space, 1, part)
end