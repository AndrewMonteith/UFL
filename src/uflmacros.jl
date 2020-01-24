export ufl_type

const EXPORT_UFL_PROPERTY_METHODS = false

field(sym::Symbol, t) = Expr(:(::), sym, t)

fields = Dict(
    :ufl_shape => (field(:ufl_shape, DimensionTuple), ()),
    :ufl_free_indices => (field(:ufl_free_indices, VarTuple{Index}), ()),
    :ufl_index_dimensions => (field(:ufl_index_dimensions, DimensionTuple), ()),
    :ufl_operands => (field(:ufl_operands, VarTuple{AbstractExpr}), ()),

    # any field marked with Any means we don't have an appropriate data type for it 
    :ufl_domain => (field(:ufl_domain, Any), ()), # This field is meant to be depreceated?
)

macro load_ufl_property_methods()
    methods = [] 
    
    for ufl_prop in keys(fields)   
        # hacky way to stop the name mangling...     
        prop_str = String(ufl_prop)
        push!(methods, esc(:($ufl_prop(x::AbstractExpr) = fields[Symbol($prop_str)][2])))
        
        if EXPORT_UFL_PROPERTY_METHODS
            push!(methods, Expr(:export, ufl_prop))
        end 
    end 

    return Expr(:block, methods...) 
end 

@load_ufl_property_methods


tag_methods = Dict(
    :inherit_indices_from_operand => (struct_name, operand_id) -> esc(quote 
        ufl_free_indices(x::$struct_name) = ufl_free_indices( ufl_operands(x)[$operand_id] )
        ufl_index_dimensions(x::$struct_name) = ufl_index_dimensions( ufl_operands(x)[$operand_id] )
    end),

    :inherit_shape_from_operand => (struct_name, operand_id) -> esc(quote 
        ufl_shape(x::$struct_name) = ufl_shape( ufl_operands(x)[$operand_id] )
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

typecode = 0

"""
    Inserts predefined fields for a struct 
    Generates accessors for the fields 
    Inserts the fields at the end of the struct body, so will be the last parameters in the new method
    If the struct contains any tags, we insert the specified additional behaviour
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
        deleteat!(struct_fields, fields_index)

        wanted_fields = []

        for wanted_field in ufl_fields.args[2].args
            field_name = Symbol(:ufl_, wanted_field)
    
            if !(field_name in keys(fields))
                error("don't recognise field $(field_name)")
            end
    
            push!(wanted_fields, fields[field_name][1])
            push!(added_methods, method(struct_name, field_name))
        end

        prepend!(struct_fields, wanted_fields)
    end 

    # if the user has specified ufl_tags field 
    # we added a method for each tag they gave
    tags_index, ufl_tags = find_field(struct_fields, :ufl_tags)
    if tags_index !== nothing 
        deleteat!(struct_fields, tags_index)

        for tag in ufl_tags.args[2].args
            if isa(tag, Symbol) 
                push!(added_methods, tag_methods[tag](struct_name))
            else
                push!(added_methods, tag_methods[tag.args[1]](struct_name, tag.args[2]))
            end
        end
    end

    tc = global typecode += 1

    push!(added_methods, esc(quote 
        ufl_typecode(x::$struct_name) = $tc
    end))

    return Expr(:block, expr, added_methods...)
end

macro attach_hash_operators(e)
    e.head === :struct || error("can only define hash operators on structs")
    struct_name = e.args[2]

    esc(quote
        $e

        Base.hash(x::$struct_name) = hash(hash_data(x))
        Base.:(==)(x::$struct_name, y::$struct_name) = hash_data(x) === hash_data(y) 
    end)
end

hash_data(n::Nothing) = nothing