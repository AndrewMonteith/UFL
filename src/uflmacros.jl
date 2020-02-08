export ufl_type, attach_hash_operators

field(sym::Symbol, t) = Expr(:(::), sym, t)

fields = Dict(
    :ufl_shape => (type=:DimensionTuple, default_val=()),
    :ufl_operands => (type=VarTuple{AbstractExpr}, default_val=()),
    :ufl_free_indices => (type=:(VarTuple{Index}), default_val=()),
    :ufl_index_dimensions => (type=:DimensionTuple, default_val=()),
)
        
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
"""
macro ufl_type(expr)
    struct_fields = expr.args[3].args

    methods_to_add, fields_to_add = [], []
    struct_name = get_struct_name(expr)
   
    fields_index, ufl_fields = find_field(struct_fields, :ufl_fields)
    if fields_index !== nothing
        deleteat!(struct_fields, fields_index)

        wanted_fields = ufl_fields.args[2].args

        for wanted_field ∈ wanted_fields
            field_name = Symbol(:ufl_, wanted_field)
    
            if !(field_name in keys(fields))
                error("don't recognise field $(field_name)")
            end
    
            push!(fields_to_add, field(field_name, fields[field_name].type))
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