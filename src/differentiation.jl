export grad 

abstract type AbstractDifferential <: Operator end

# abstract type CompoundDerivative <: AbstractDifferential end 

@ufl_type struct Grad <: AbstractDifferential
    ufl_fields = (operands,shape)
    ufl_tags = (num_ops=1,)

    dim::Dimension

    function Grad(f)
        is_cellwise_constant(f) && return Zero(tuple(ufl_shape(f)..., geometric_dimension(f)), ufl_free_indices(f), ufl_index_dimensions(f))

        gdim = find_geometric_dimension(f)
        new(@sig((f,)), tuple(ufl_shape(f)..., gdim), gdim)
    end
end
Base.show(io::IO, g::Grad) = print(io, "grad($(g.ufl_operands[1]))")

grad(f) = (Grad ∘ as_ufl)(f)