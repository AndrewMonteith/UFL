export Constant

constant_id = 0

@ufl_type struct Constant 
    ufl_fields = (domain, shape)

    id::Int 

    function Constant(domain, shape::DimensionTuple=())
        id = constant_id
        global constant_id = constant_id + 1

        new(id, domain, shape)
    end
end

is_cellwise_constant(::Constant) = true 
ufl_domains(c::Constant) = (ufl_domain(c),)
Base.repr(io::IO, c::Constant) = "Constant(domain $(c.ufl_domain) shape $(c.ufl_shape) count $(c.id)"
Base.show(io::IO, c::Constant) = show(io, "c_$(c.id)")