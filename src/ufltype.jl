export ufl_type

macro ufl_type()
    fields = []
    add_field! = (sym::Symbol, t)->append!(fields, Expr(:(::), esc(sym), t))

    add_field!(:ufl_free_indicies, DimensionTuple)
    add_field!(:ufl_shape, DimensionTuple)
    add_field!(:ufl_free_indicies, VarTuple{FreeIndicies})
    add_field!(:ufl_operands, VarTuple{AbstractExpr})

    esc(quote $(fields...) end)
end