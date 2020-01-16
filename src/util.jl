function get_opt_type(type)
    types = if isa(type, Symbol)
        [type]
    else
        union_types = type.args

        union_types[2:length(union_types)]
    end

    :(Union{$(types...), Nothing})
end

macro opt_t(type)
    esc(get_opt_type(type))
end

macro opt(e)
    # Takes a variable declaration and x::T and transforms it to x::Union{T, Nothing}=nothing
    # If e is a x::Union{T1, ..., TN} it becomes x::Union{T1, ..., TN, Nothing}=nothing

    sym, types = e.args[1], get_opt_type(e.args[2])

    ex = :($sym::$types=nothing)
    ex.head = :kw
    return esc(ex)
end