export Constant

constant_id = 0

@ufl_type struct Constant 
    ufl_fields = (domain, shape)

    id::Int 

    function Constant(shape::DimensionTuple=())
        id = constant_id
        global constant_id = constant_id + 1

        new(id, shape)
    end
end

is_cellwise_constant(::Constant) = true 
ufl_domains(c::Constant) = (ufl_domain(c),)
pretty_print(c::Constant) = "Constant(domain $(c.ufl_domain) shape $(c.ufl_shape) count $(c.id)"
mathsy_print(c::Constant) = "c_$(c.id)"