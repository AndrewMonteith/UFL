export ufl_type

field(sym::Symbol, t) = Expr(:(::), sym, t)

fields = Dict(
    :ufl_shape => field(:ufl_shape, DimensionTuple),
    :ufl_free_indices => field(:ufl_free_indices, VarTuple{Index}),
    :ufl_index_dimensions => field(:ufl_index_dimensions, DimensionTuple),
    :ufl_operands => field(:ufl_operands, VarTuple{AbstractExpr}),
    :ufl_domain => field(:ufl_domain, Any) 
)

tag_methods = Dict(
    :inherit_indices_from_operand => (struct_name, operand_id) -> esc(quote 
        ufl_free_indices(x::$struct_name) = ufl_free_indices(ufl_operands(x)[operand_id])
        ufl_index_dimensions(x::$struct_name) = ufl_index_dimensions(ufl_operands(x)[operand_id])
    end)
)

method(name::Symbol, param::Symbol) = esc(quote $param(x::$name) = x.$param end)
get_struct_name(def::Expr) = isa(def.args[2], Symbol) ? def.args[2] : def.args[2].args[1]

function find_field(fields::AbstractArray{Any}, field_name)
    for (i, field) in enumerate(fields)
        if isa(field, Expr) && field.args[1] == field_name 
            return i, field 
        end
    end
    
    return (nothing, nothing)
end

"""
    Inserts predefined fields for a struct 
    Generates accessors for the fields 
    Inserts the fields at the end of the struct body, so will be the last parameters in the new method
"""
macro ufl_type(expr)
    struct_fields = expr.args[3].args
    added_methods = []

    struct_name = get_struct_name(expr)
    
    # Has the user provided ufl_fields tag in struct body 
    # if so we generte all the properties they requested as members on the struct 
    # and generate an accessor for the property that just returns the member 
    fields_index, ufl_fields = find_field(struct_fields, :ufl_fields)
    if fields_index !== nothing
        for wanted_field in ufl_fields.args[2].args
            field_name = Symbol(:ufl_, wanted_field)
    
            if !(field_name in keys(fields))
                error("don't recognise field $(field_name)")
            end
    
            push!(struct_fields, fields[field_name])
            push!(added_methods, method(struct_name, field_name))
        end

        deleteat!(struct_fields, fields_index)
    end 

    # if the user has specified ufl_tags field 
    # we added a method for each tag they gave
    tags_index, ufl_tags = find_field(struct_fields, :ufl_tags)
    if tags_index !== nothing 
        for tag in ufl_tags.args[2].args
            if isa(tag, Symbol) 
                push!(added_methods, tag_methods[tag](struct_name))
            else
                push!(added_methods, tag_methods[tag.args[1]](struct_name, tag.args[2]))
            end
        end

        deleteat!(struct_fields, tags_index)
    end

    return Expr(:block, expr, added_methods...)
end
