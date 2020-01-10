export UflConstant, Constant

constant_id = 0

@ufl_type struct UflConstant 
    """ Mirrors the Constant Literal class as defined in UFL """
    ufl_fields = (domain, shape)

    id::Int 

    function Constant(domain, shape::DimensionTuple=())
        id = constant_id
        global constant_id = constant_id + 1

        new(domain, shape, id)
    end
end

is_cellwise_constant(::UflConstant) = true 
ufl_domains(c::UflConstant) = (ufl_domain(c),)
Base.repr(io::IO, c::UflConstant) = "Constant(domain $(c.ufl_domain) shape $(c.ufl_shape) count $(c.id)"
Base.show(io::IO, c::UflConstant) = show(io, "c_$(c.id)")



struct Constant  
    """ Constant literal as defined in the firedrake proejct """
    
    
end



























