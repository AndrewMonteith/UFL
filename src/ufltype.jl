export ufl_type

"""
    This macro may not be needed to be honest
"""

field(sym::Symbol, t) = Expr(:(::), sym, t)

fields = Dict(
    :ufl_shape => field(:ufl_shape, DimensionTuple),
    :ufl_free_indices => field(:ufl_free_indices, VarTuple{Index}),
    :ufl_index_dimensions => field(:ufl_index_dimensions, DimensionTuple),
    :ufl_operands => field(:ufl_operands, VarTuple{AbstractExpr}),
    :ufl_domain => field(:ufl_domain, Any) 
)

macro ufl_type(symbols...)
    struct_symbol = symbols[1]
    field_exprs = []

    print(symbols)

    for sym in symbols[2:length(symbols)-1] 
        println("Looking for $(sym)")
        println("Inside:", keys(fields))
        if !(sym in keys(fields))
            error("$(sym) does not exist in fields dict")
        end

        # print("Loading symbol:", )
        # eval(quote 
        #     $sym(x::$struct_symbol) = x.$sym 
        # end)

        push!(field_exprs, fields[field])
    end 

    esc(quote $(field_exprs...) end)
end
