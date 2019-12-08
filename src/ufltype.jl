export ufl_type

"""
    Master macro for creating UFL types
    Currently just puts some fields on the struct 
    In future will make it do some checks in the ctor such as:
        - len(ufl_free_indices) === len(ufl_index_dimensions)
"""
macro ufl_type(has_operands::Bool=false)
    fields = Vector{Expr}()
    add_field! = (sym::Symbol, t) -> push!(fields, Expr(:(::), sym, t))

    add_field!(:ufl_shape, DimensionTuple)
    add_field!(:ufl_free_indices, VarTuple{Index})
    add_field!(:ufl_index_dimensions, DimensionTuple)

    if has_operands
        add_field!(:ufl_operands, VarTuple{AbstractExpr})
    end

    esc(quote $(fields...) end)
end