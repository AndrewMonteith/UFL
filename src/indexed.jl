export Indexed 

struct Indexed <: Operator 
    ufl_operands::Tuple{AbstractExpr, MultiIndex}

    function Indexed(expression::AbstractExpr, m::MultiIndex)
        new((expression, m))
    end

end