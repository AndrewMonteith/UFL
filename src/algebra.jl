export Sum 

@ufl_type struct Sum <: Operator 
    ufl_fields = (operands,)
    ufl_tags = (inherit_shape_from_operand=1, inherit_indices_from_operand=1)
    
    function Sum(a::AbstractExpr, b::AbstractExpr)
        shape = ufl_shape(a)
        if shape !== ufl_shape(b) 
            error("Cannot add expressions of different shape")
        end 

        free_indicies = ufl_free_indices(a)
        if free_indicies != ufl_free_indices(b) 
            error("Cannot add indices with different type")
        end 

        index_dimensions = ufl_index_dimensions(a)
        if index_dimensions != ufl_index_dimensions(b)
            error("Cannot add expressions with different index dimensions")
        end 

        if a isa Zero 
            return b 
        elseif b isa Zero 
            return a 
        end 

        a_is_scalar, b_is_scalar = a isa ScalarValue, b isa ScalarValue 

        if a_is_scalar && b_is_scalar
            return ScalarValue(a.val + b.val)
        elseif b_is_scalar
            a, b = b, a
        else
            # Do some kind of convoluted sorting algorithm which I hope I don't 
            # have to implement 
        end 

        new((a, b))
    end
end 


Base.:+(e1::AbstractExpr, e2::AbstractExpr) = Sum(e1, e2)