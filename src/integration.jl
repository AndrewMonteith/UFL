export Measure, Form

measure_name_to_integral_type = Dict("dx" => "cell", "ds" => "exterior_facet")
integral_type_to_measure_name = Dict(j => i for (i, j) ∈ measure_name_to_integral_type)

const Subdomain = Union{String, Int}

struct Measure 
    domain::@opt_t(Mesh)
    integral_type::String 
    metadata::Dict{Any, Any} 
    subdomain_id::Subdomain
    subdomain_data 
    
    function Measure(measure_name::String)
        integral_type = measure_name_to_integral_type[measure_name]
        new(nothing, integral_type, Dict(), "everywhere", nothing)
    end

    function Measure(domain::@opt_t(Mesh), integral_type::String, metadata::Dict{Any, Any}, subdomain_id::Subdomain, subdomain_data)
        new(domain, integral_type, metadata, subdomain_id, subdomain_data)
    end 
end 
Base.hash(m::Measure) = hash((m.integral_type, m.subdomain_id, m.subdomain_data, m.metadata, m.domain))
Base.show(io::IO, m::Measure) = print(io, integral_type_to_measure_name[m.integral_type])

function reconstruct(x::T; kwargs...) where T
    fields = fieldnames(T)
    new_members = []

    for field ∈ fields
        push!(new_members,  get(kwargs, field, getfield(x, field)))
    end

    T(new_members...)
end 



function (m::Measure)(@opt(subdomain_id::Union{Int, String, Mesh}), @opt(domain::Mesh))
    # if no args are provided
    subdomain_id === nothing && reconstruct(m; subdomain_id="everywhere") 

    if subdomain_id isa Mesh 
        subdomain_id, domain = "everywhere", subdomain_id
    end 

    reconstruct(m; subdomain_id=subdomain_id, domain=domain)
end

function join_domains(domains::Vector{Mesh})::VarTuple{Mesh}
    isempty(domains) && return () 

    unique!(geometric_dimension, domains)

    length(domains) > 1 && error("Multiple domains with different geometric dimensions")

    tuple(domains...)
end

function extract_domains(expr::AbstractExpr)::VarTuple{Mesh}
    domains::Vector{Mesh} = [] 

    @UFL.unique_pre_traversal for op ∈ expr
        d = ufl_domain(op) 
        if d !== nothing 
            push!(domains, d)
        end
    end

    join_domains(domains)
end 

function Base.:*(integrand::AbstractExpr, m::Measure)
    !is_true_scalar(integrand) && error("Can only integrate over scalar expressions.")

    # In Future: Might need to support subdomain_id as a tuple?
    domain = m.domain 
    if domain === nothing 
        domains = extract_domains(integrand)
        domain = if length(domains) === 1
            domains[1]
        elseif length(domains) === 0
            error("Integral is missing an integration domain")
        else 
            error("Multiple domains are present")
        end 
    end 


    integral = Integral(integrand, m.integral_type, domain, m.subdomain_id)

    Form((integral,))
end

Base.:*(expr::Real, m::Measure) = as_ufl(expr)*m


struct Integral
    integrand::AbstractExpr
    integral_type::String
    domain::Mesh 
    subdomain_id::Subdomain
    
    function Integral(integrand::AbstractExpr, integral_type::String, domain::Mesh, subdomain_id::Subdomain)
        new(integrand, integral_type, domain, subdomain_id)
    end
end

Base.:-(i::Integral) = reconstruct(i; integrand=-i.integrand)

function Base.show(io::IO, i::Integral)
    measure_name = integral_type_to_measure_name[i.integral_type]
    
    print(io, "{ $(i.integrand) } * $(measure_name)($(i.domain)[nothing])")
end

struct Form  
    integrals::VarTuple{Integral}
    
    function Form(integrals::VarTuple{Integral})
        new(integrals)
    end 
end 

Base.:+(f::Form, f2::Form) = Form(tuple(f.integrals..., f2.integrals...))
Base.:+(f::Form, r::Real) = r === 0 ? f : error("cannot add non-zero real to a form")
Base.:+(r::Real, f::Form) = f+r
function Base.:+(f::Form, zero::Zero)
    if ((isempty ∘ ufl_shape)(zero) && (isempty ∘ ufl_free_indices))
        zero
    else
        error("cannot add non-scalar zero to a form")
    end
end
Base.:+(zero::Zero, f::Form) = f+zero
Base.:-(f::Form, f2::Form) = f + -f2
Base.:-(f::Form, x) = f + (-x)
Base.:-(x, f::Form) = (-x) + f
Base.:-(f::Form) = Form(tuple((-integral for integral ∈ f.integrals)...))

function Base.:*(f::Form, s::AbstractExpr)
    if is_scalar_constant_expression(s)
        Form([s*integral for integral ∈ f.integrals])
    else
        error("still need to implement action form.py:326")
    end
end 