export ufl_type, ufl_shape, ufl_free_indices, ufl_index_dimensions

using OrderedCollections: LittleDict

field(sym::Symbol, t) = Expr(:(::), sym, t)
method(name::Symbol, param::Symbol) = esc(quote $param(x::$name) = x.$param end)

required_fields = LittleDict(
    :ufl_shape => (type=:DimensionTuple, default_val=()),
    :ufl_free_indices => (type=:(VarTuple{Index}), default_val=()),
    :ufl_index_dimensions => (type=:DimensionTuple, default_val=()),
)

optional_fields = Dict(
    :ufl_operands => (type=VarTuple{AbstractExpr}, default_val=())
)

tag_methods = Dict(
    :inherit_indices_from_operand => (struct_name, operand_id) -> esc(quote 
        ufl_free_indices(x::$struct_name) = ufl_free_indices( ufl_operands(x)[$operand_id] )
        ufl_index_dimensions(x::$struct_name) = ufl_index_dimensions( ufl_operands(x)[$operand_id] )
    end),

    :inherit_shape_from_operand => (struct_name, operand_id) -> esc(quote 
        ufl_shape(x::$struct_name) = ufl_shape( ufl_operands(x)[$operand_id] )
    end)
)

for (required_field, def) ∈ required_fields 
    @eval begin 
        $required_field(x::AbstractExpr) = x.$required_field 
        export $required_field
    end 
end

for (optional_field, def) ∈ optional_fields 
    @eval begin 
        $optional_field(x::AbstractExpr) = $(def.default_val)
    end
end
        
function find_field(fields::AbstractArray{Any}, field_name::Symbol)
    for (i, field) in enumerate(fields)
        if isa(field, Expr) && field.args[1] == field_name 
            return i, field 
        end
    end
    
    return (nothing, nothing)
end

typecode = 0

function get_struct_name(def::Expr)
    if def.args[2] isa Symbol  # Struct that does not subtype abstract type
        def.args[2]
    elseif def.args[2].args[1] isa Symbol # Subtyped struct with no parametric types 
        def.args[2].args[1]
    else
        def.args[2].args[1].args[1]
    end
end

"""
    Inserts predefined fields for a struct 
    Generates accessors for the fields 
    Inserts the fields at the end of the struct body, so will be the last parameters in the new method
    If the struct contains any tags, we insert the specified additional behaviour
"""
macro ufl_type(expr)
    struct_fields = expr.args[3].args

    methods_to_add, fields_to_add = [], []
    struct_name = get_struct_name(expr)
   
    # Any struct that uses this macro will have any field from required_fields injected into it 
    # In the order specified in the dictionary
    for (name, details) ∈ required_fields 
        push!(fields_to_add, field(name, details.type))
    end

    # If the user has provided a ufl_fields tag in the struct body 
    # Any property mentioned in there must be in the optional_fields dictionary 
    # We then generate an accessor that reads from the field of the struct rather than the 
    # the default value 
    fields_index, ufl_fields = find_field(struct_fields, :ufl_fields)
    if fields_index !== nothing
        deleteat!(struct_fields, fields_index)

        for wanted_field ∈ ufl_fields.args[2].args
            field_name = Symbol(:ufl_, wanted_field)
    
            if !(field_name in keys(optional_fields))
                error("don't recognise field $(field_name)")
            end
    
            push!(fields_to_add, field(field_name, optional_fields[field_name].type))

            if field_name !== :ufl_operands
                push!(!methods_to_add, esc(:( $field_name(x::$struct_name)=x.$field_name )))
            end
        end
    end 

    # if the user has specified ufl_tags field 
    # we added a method for each tag they gave
    tags_index, ufl_tags = find_field(struct_fields, :ufl_tags)
    if tags_index !== nothing 
        deleteat!(struct_fields, tags_index)

        for tag in ufl_tags.args[2].args
            if isa(tag, Symbol) 
                push!(methods_to_add, tag_methods[tag](struct_name))
            else
                push!(methods_to_add, tag_methods[tag.args[1]](struct_name, tag.args[2]))
            end
        end
    end

    prepend!(struct_fields, fields_to_add)
    
    tc = global typecode += 1
    push!(methods_to_add, esc(quote ufl_typecode(x::$struct_name) = $tc end))

    return Expr(:block, expr, methods_to_add...)
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