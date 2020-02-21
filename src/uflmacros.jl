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

function inject_hash_behaviour(expr, typecode)
    # Assumes last line for the inner constructor is a new

    functions = findall(expr -> expr isa Expr && expr.head === :function, expr.args[3].args)
    for func_i in functions 
        inner_ctor = expr.args[3].args[func_i]

        new_call = inner_ctor.args[2].args[end].args
        
        sig_param_is = findall(param -> param isa Expr && param.head === :macrocall && param.args[1] === Symbol("@sig"), new_call)
        sig_exprs = []
        for sig_param_i ∈ sig_param_is
            # sig_param_i will be the index of parameter wrapped in a @sig 
            # we unwrap it and mark that expression needs to be included in the hash
            new_call[sig_param_i] = new_call[sig_param_i].args[3]
            push!(sig_exprs, new_call[sig_param_i])
        end

        hash_expr = :( $(esc(:compute_hash))($typecode, $(sig_exprs...)) )

        insert!(new_call, 2, hash_expr)
    end
end

# The below function was an attempt at an optimization for inner constructors of the form:
#      new(compute_hash(..., e_1, e_2, e_3), e_1, e_2, e_3)  
#
# To be optimised to:
#      h1 = e_1 
#      h2 = e_2 
#      h3 = e_3 
#      new(compute_hash(..., h1, h2, h3), h1, h2, h3) 
#
#  This theoretically would mean that the expressions aren't recomputed twice and you know what
#  LLVM does this for you already. So all this code is basically useless for now. If ever e_k as 
#  some function with a sideeffect perhaps LLVM would not perform this optimisation and this function 
#  might become useful
#
# function inject_hash_behaviour(expr, typecode)
#     # println("Loading expression: ", expr)
#     h_id = 0
#     new_hash_identifier() = Symbol(:hash, h_id += 1)
#
#     inner_ctor_is = findall(expr -> expr isa Expr && expr.head === :function, expr.args[3].args) 
#     for inner_ctor_i in inner_ctor_is 
#         inner_ctor = expr.args[3].args[inner_ctor_i]
#
#         new_call = inner_ctor.args[2].args[end].args
#         
#         sig_param_is = findall(param -> param isa Expr && param.head === :macrocall && param.args[1] === Symbol("@sig"), new_call)
#
#         hash_identifiers::Vector{Symbol} = []
#         hash_assignments::Vector{Expr} = []
#         for sig_param_i ∈ sig_param_is
#             # sig_param_i will be the index of parameter wrapped in a @sig 
#             # we unwrap it and mark that expression needs to be included in the hash
#             new_call[sig_param_i] = new_call[sig_param_i].args[3]
#
#             hash_id = new_hash_identifier()
#             push!(hash_identifiers, hash_id)            
#             push!(hash_assignments, :($(esc(hash_id)) = $(new_call[sig_param_i])))
#         end
#
#         hash_expr = (:(compute_hash($typecode, $(hash_identifiers...))))
#
#         if length(hash_assignments) === 1 
#             insert!(inner_ctor.args[2].args, length(inner_ctor.args[2].args)-1, hash_assignments...)
#         elseif length(hash_assignments) > 1 
#             insert!(inner_ctor.args[2].args, length(inner_ctor.args[2].args)-1, hash_assignments)
#         end
#
#         insert!(new_call, 2, esc(hash_expr))
#     end
# end

function insert_reconstructor(expr, struct_name, struct_fields)
    member_variables = filter(line -> line isa Expr && line.head === :(::), struct_fields)

    !any(var -> var.args[1] === :ufl_operands, member_variables) && return

    rector_new_args = map(var -> var.args[1] === :ufl_operands ? :operands : :(x.$(var.args[1])), member_variables)

    recontructor = esc(:(
        function $struct_name(x::$struct_name, operands::VarTuple{UFL.AbstractExpr}) 
            new($(rector_new_args...))
        end
    ))

    push!(struct_fields, recontructor)
end

"""
    Inserts predefined fields for a struct
    Parses ufl_fields tags into metadata
    Insert copy constructor
    Injects hash logic
"""
macro ufl_type(expr)
    struct_fields = expr.args[3].args

    methods_to_add, fields_to_add = [], []
    struct_name = get_struct_name(expr)

    prepend!(fields_to_add, (field(:ufl_hash_code, :UInt32),))

    metadata = Dict{Symbol, Int}()
    tags_index, ufl_tags = find_field(struct_fields, :ufl_tags)
    if tags_index !== nothing 
        deleteat!(struct_fields, tags_index)

        wanted_tags = ufl_tags.args[2].args
        for wanted_tag ∈ wanted_tags 
            metadata[wanted_tag.args[1]] = wanted_tag.args[2]
        end
    end

    # Expand ufl_fields variable into respective struct members
    fields_index, ufl_fields = find_field(struct_fields, :ufl_fields)
    if fields_index !== nothing
        deleteat!(struct_fields, fields_index)

        wanted_fields = ufl_fields.args[2].args

        for wanted_field ∈ wanted_fields
            field_name = Symbol(:ufl_, wanted_field)
    
            field_name ∉ keys(fields) && error("don't recognise field $(field_name)")

            if field_name === :ufl_operands 
                num_operands = get!(metadata, :num_ops, -1)

                var_type = num_operands === -1 ? fields[field_name].type : :(NTuple{$(num_operands), AbstractExpr})
                push!(fields_to_add, field(field_name, var_type))
            else
                push!(fields_to_add, field(field_name, fields[field_name].type))
            end
        end
    end 

    # Insert Copy Constructor 
    tc = global typecode += 1
    push!(methods_to_add, esc(quote ufl_typecode(x::$struct_name) = $tc end))
    
    prepend!(struct_fields, fields_to_add)

    insert_reconstructor(expr, struct_name, struct_fields)
    inject_hash_behaviour(expr, tc)

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