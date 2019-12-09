field(sym::Symbol, t) = Expr(:(::), sym, t)
method(name::Symbol, param::Symbol) = esc(quote 
    $param(x::$name) = x.$param 
end)

macro foobar(expr)
    struct_name = expr.args[2]
    struct_fields = expr.args[3].args
    struct_methods = []
    
    ufl_fields = struct_fields[2]

    if ufl_fields.args[1] != :ufl_fields 
        error("malformed UFL type")
    end
   
    for wanted_field in ufl_fields.args[2].args
        field_name = Symbol(:ufl_, wanted_field)

        push!(struct_fields, field(field_name, Int))
        push!(struct_methods, method(struct_name, field_name))
    end

    deleteat!(struct_fields, 2)

    new_expr = Expr(:block, expr, struct_methods...)

    return new_expr
end

@foobar struct X 
    ufl_fields = (shape,)
end

x = X(1)

ufl_shape(x)