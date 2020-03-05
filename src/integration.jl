export Measure

measure_name_to_integral_type = Dict("dx" => "cell", "ds" => "exterior_facet")
integral_type_to_measure_name = Dict(j => i for (i, j) ∈ measure_name_to_integral_type)

struct Measure 
    domain 
    integral_type::String 
    metadata::Dict{Any, Any} 
    subdomain_id::Union{String, Int}
    subdomain_data 
    
    function Measure(measure_name::String)
        integral_type = measure_name_to_integral_type[measure_name]
        new(nothing, integral_type, Dict(), "everywhere", nothing)
    end
end 
Base.hash(m::Measure) = hash((m.integral_type, m.subdomain_id, m.subdomain_data, m.metadata, m.domain))
Base.show(io::IO, m::Measure) = print(io, integral_type_to_measure_name[m.integral_type])

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
    subdomain_id::String 
    
    function Integral(integrand::AbstractExpr, integral_type::String, domain::Mesh, subdomain_id::String)
        new(integrand, integral_type, domain, subdomain_id)
    end
end

function Base.show(io::IO, i::Integral)
    measure_name = integral_type_to_measure_name[i.integral_type]

    print(io, "{ $(i.integrand) } * $(measure_name)($(i.domain)[$(i.subdomain_data)])")
end


struct Form  
    integrals::VarTuple{Integral}

    function Form(integrals::VarTuple{Integral})
        new(integrals)
    end 
end 